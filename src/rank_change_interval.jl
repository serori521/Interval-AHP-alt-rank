module RankChangeInterval
# 順序変化点を求めるモジュール
using Statistics, LinearAlgebra
export find_rank_change_points, analyze_all_alternatives
include("calc_IPW.jl")

# 指定した t における代替案 o_o と o_q 間のリグレット値を返す関数
function regret_value(o_o::Int, o_q::Int, t::Float64, matrix, methodW)
    if o_o == o_q
        return 0.0
    end
    
    Y_L = methodW.L .* t
    Y_R = methodW.R .* t
    
    # 一時的に regret を計算
    temp_matrix = deepcopy(matrix)
    _, _ = calc_regret(temp_matrix, Y_L, Y_R)
    
    # 行列内のリグレット値を取得
    i, j = minmax(o_o, o_q)
    if i == o_o
        return temp_matrix[i, j].regret
    else
        return -temp_matrix[i, j].regret  # 符号を反転
    end
end

# 2点評価から直線： r(t) = A*t + B を返す（t0,t1 における評価値から求める）
function linear_regret_parameters(o_o::Int, q::Int, t0::Float64, t1::Float64, matrix, methodW)
    r0 = regret_value(o_o, q, t0, matrix, methodW)
    r1 = regret_value(o_o, q, t1, matrix, methodW)
    
    A = (r1 - r0) / (t1 - t0)
    B = r0 - A*t0
    return A, B
end

# 2候補（p と q）のリグレット関数の直線同士の交点 t を求める
function crossing_point(o_o::Int, p::Int, q::Int, t0::Float64, t1::Float64, matrix, methodW)
    A_p, B_p = linear_regret_parameters(o_o, p, t0, t1, matrix, methodW)
    A_q, B_q = linear_regret_parameters(o_o, q, t0, t1, matrix, methodW)
    
    if isapprox(A_p, A_q; atol=1e-8)
        # 傾きがほぼ同じ場合は、中間値を返す
        return (t0 + t1)/2
    else
        t_cross = (B_q - B_p) / (A_p - A_q)
        # [t0, t1] 内に収める
        t_cross = max(t0, min(t_cross, t1))
        return t_cross
    end
end

# 指定の対象代替案 o_o に対して、区間 [t0, t1] 内で順序が変わる t の境界を求める
function find_rank_change_points(o_o::Int, t0::Float64, t1::Float64, 
                                matrix, methodW; 
                                max_iter::Int=20, tol::Float64=1e-6)
    n = size(matrix, 1)
    candidates = [q for q in 1:n if q != o_o]
    
    # 初期化: t0 での最大リグレットを与える代替案を見つける
    regrets_at_t0 = [regret_value(o_o, q, t0, matrix, methodW) for q in candidates]
    p_idx = argmax(regrets_at_t0)
    p = candidates[p_idx]
    
    # 2番目に大きいリグレットを与える代替案
    temp_regrets = copy(regrets_at_t0)
    temp_regrets[p_idx] = -Inf
    u_idx = argmax(temp_regrets)
    u = candidates[u_idx]
    
    # t1 でのリグレットを計算
    regrets_at_t1 = [regret_value(o_o, q, t1, matrix, methodW) for q in candidates]
    
    # t1 で p より大きいリグレットを与える候補を見つける
    O = [candidates[i] for i in 1:length(candidates) 
         if regrets_at_t1[i] > regrets_at_t1[p_idx]]
    
    # 各候補との交点を計算
    t_cross = Dict{Int, Float64}()
    for q in O
        t_cross[q] = crossing_point(o_o, p, q, t0, t1, matrix, methodW)
    end
    
    # 不要な候補を除外
    for q in copy(O)
        for s in O
            if s != q && t_cross[s] < t_cross[q] && 
               regret_value(o_o, s, t1, matrix, methodW) > regret_value(o_o, q, t1, matrix, methodW)
                filter!(x -> x != q, O)
                break
            end
        end
    end
    
    # 結果を格納する配列
    change_points = Float64[]
    max_regrets = Float64[]
    max_regret_alts = Int[]
    
    # 反復処理
    iter = 0
    current_t0 = t0
    current_p = p
    
    while iter < max_iter && current_t0 < t1 - tol
        # 次の交点を見つける
        next_t = t1
        next_q = 0
        
        for q in O
            if t_cross[q] > current_t0 && t_cross[q] < next_t
                next_t = t_cross[q]
                next_q = q
            end
        end
        
        if next_q == 0 || isapprox(next_t, current_t0, atol=tol)
            # 次の交点がない、または現在の位置と同じ場合は終了
            break
        end
        
        # 交点を記録
        push!(change_points, next_t)
        
        # 交点でのリグレットを計算
        regret_at_cross = regret_value(o_o, next_q, next_t, matrix, methodW)
        push!(max_regrets, regret_at_cross)
        push!(max_regret_alts, next_q)
        
        # 更新
        current_t0 = next_t
        current_p = next_q
        
        iter += 1
    end
    
    # 最終的な候補集合と最大・最小リグレット
    O_hat = union(O, [p])
    r_M = maximum([regret_value(o_o, q, t1, matrix, methodW) for q in O_hat])
    r_m = minimum([regret_value(o_o, q, t0, matrix, methodW) for q in O_hat])
    
    return (
        change_points = change_points,
        max_regrets = max_regrets,
        max_regret_alts = max_regret_alts,
        O_hat = O_hat,
        r_M = r_M,
        r_m = r_m
    )
end

# すべての代替案について順位変化点を分析
function analyze_all_alternatives(utility::Matrix{Float64}, methodW, t_range::Tuple{Float64, Float64})
    n = size(utility, 1)
    matrix = create_minimax_R_Matrix(utility)
    
    # 各代替案の順位変化点を格納
    all_results = Dict{Int, Any}()
    
    for o_o in 1:n
        result = find_rank_change_points(o_o, t_range[1], t_range[2], matrix, methodW)
        all_results[o_o] = result
    end
    
    # 全ての変化点をソート
    all_change_points = Float64[]
    for o_o in 1:n
        append!(all_change_points, all_results[o_o].change_points)
    end
    unique!(sort!(all_change_points))
    
    # 各変化点での順位を計算
    rankings = Dict{Float64, Vector{Int}}()
    
    # 初期点と終点を追加
    t_points = [t_range[1]; all_change_points; t_range[2]]
    
    for t in t_points
        # 各代替案の最大リグレットを計算
        max_regrets = zeros(n)
        for o_o in 1:n
            # t での最大リグレットを計算
            max_r = -Inf
            for q in 1:n
                if q != o_o
                    r = regret_value(o_o, q, t, matrix, methodW)
                    max_r = max(max_r, r)
                end
            end
            max_regrets[o_o] = max_r
        end
        
        # リグレットの小さい順に並べる（順位付け）
        rankings[t] = sortperm(max_regrets)
    end
    
    return (
        all_results = all_results,
        change_points = all_change_points,
        rankings = rankings
    )
end

end # module
