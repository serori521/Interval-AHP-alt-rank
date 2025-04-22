module RankChangeInterval
using Statistics, LinearAlgebra
export find_rank_change_points_from_tR, analyze_all_alternatives_with_avail_space
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

# アベイルスペースから次のt'を計算する関数
function compute_next_t(current_t::Float64, matrix, methodW)
    Y_L = methodW.L .* current_t
    Y_R = methodW.R .* current_t
    
    temp_matrix = deepcopy(matrix)
    _, _ = calc_regret(temp_matrix, Y_L, Y_R)
    
    # 非対角要素の最小アベイルスペースを見つける
    min_availspace = findmin([temp_matrix[i,j].Avail_space 
                             for i in 1:size(matrix, 1), j in 1:size(matrix, 1) 
                             if i != j])[1]
    
    # rパラメータを計算
    r = 1 / (1 + min_availspace)
    
    # 次のtを計算
    next_t = r * current_t
    
    return next_t, min_availspace
end

# t内でのリグレット関数の線形パラメータ（傾きと切片）を計算
function linear_regret_parameters(o_o::Int, q::Int, t1::Float64, t2::Float64, matrix, methodW)
    r1 = regret_value(o_o, q, t1, matrix, methodW)
    r2 = regret_value(o_o, q, t2, matrix, methodW)
    
    # 線形関数 r(t) = A*t + B のパラメータを計算
    A = (r2 - r1) / (t2 - t1)
    B = r1 - A*t1
    
    return A, B
end

# 2つのリグレット関数の交点を計算
function compute_crossing_point(o_o::Int, p::Int, q::Int, t_upper::Float64, t_lower::Float64, matrix, methodW)
    A_p, B_p = linear_regret_parameters(o_o, p, t_upper, t_lower, matrix, methodW)
    A_q, B_q = linear_regret_parameters(o_o, q, t_upper, t_lower, matrix, methodW)
    
    if isapprox(A_p, A_q; atol=1e-8)
        # 傾きがほぼ同じ場合は交点なしとみなす
        return nothing
    end
    
    # 交点のt値を計算
    t_cross = (B_q - B_p) / (A_p - A_q)
    
    # t_upperとt_lowerの間にあるか確認
    if t_cross >= t_lower && t_cross <= t_upper
        return t_cross
    else
        return nothing
    end
end

# 指定した区間[t_upper, t_lower]内での順位変化点を検出
function find_rank_change_in_interval(o_o::Int, t_upper::Float64, t_lower::Float64, matrix, methodW)
    n = size(matrix, 1)
    candidates = [q for q in 1:n if q != o_o]
    
    # t_upperでの最大リグレットを与える代替案を特定
    regrets_at_upper = [regret_value(o_o, q, t_upper, matrix, methodW) for q in candidates]
    p_idx = argmax(regrets_at_upper)
    p = candidates[p_idx]
    
    # t_lowerでp以外に最大リグレットを与える可能性がある代替案を特定
    regrets_at_lower = [regret_value(o_o, q, t_lower, matrix, methodW) for q in candidates]
    potential_winners = [candidates[i] for i in 1:length(candidates) 
                        if regrets_at_lower[i] >= regrets_at_lower[p_idx]]
    
    # 交点と順位変化点を記録
    change_points = Float64[]
    max_regret_alts = Int[]
    
    for q in potential_winners
        if q == p
            continue
        end
        
        t_cross = compute_crossing_point(o_o, p, q, t_upper, t_lower, matrix, methodW)
        if t_cross !== nothing
            # この交点でのリグレットを計算
            r_cross = regret_value(o_o, q, t_cross, matrix, methodW)
            
            # 本当にこの時点で最大リグレットを与えるか確認
            is_max = true
            for other in candidates
                if other != q && regret_value(o_o, other, t_cross, matrix, methodW) > r_cross
                    is_max = false
                    break
                end
            end
            
            if is_max
                push!(change_points, t_cross)
                push!(max_regret_alts, q)
            end
        end
    end
    
    return change_points, max_regret_alts
end

# t^Rから開始して、アベイルスペースに基づいて順位変化点を検出
function find_rank_change_points_from_tR(o_o::Int, t_L::Float64, t_R::Float64, matrix, methodW)
    # 初期化
    change_points = Float64[]
    max_regret_alts = Int[]
    intervals = Tuple{Float64, Float64}[]
    
    # 現在のt値を初期化
    current_t = t_R
    
    while current_t > t_L
        # 次のt値を計算
        next_t, _ = compute_next_t(current_t, matrix, methodW)
        
        # 下限t_Lを下回らないようにする
        next_t = max(next_t, t_L)
        
        # 区間を記録
        push!(intervals, (next_t, current_t))
        
        # 区間内の順位変化点を検出
        interval_changes, interval_alts = find_rank_change_in_interval(o_o, current_t, next_t, matrix, methodW)
        
        # 全体の結果に追加
        append!(change_points, interval_changes)
        append!(max_regret_alts, interval_alts)
        
        # 次の区間へ
        current_t = next_t
        
        # t_Lに到達したら終了
        if isapprox(current_t, t_L, atol=1e-8)
            break
        end
    end
    
    # 各区間での最大・最小リグレット値を計算
    r_M = regret_value(o_o, argmax([regret_value(o_o, q, t_L, matrix, methodW) 
                                  for q in 1:size(matrix, 1) if q != o_o]), 
                       t_L, matrix, methodW)
    
    r_m = regret_value(o_o, argmax([regret_value(o_o, q, t_R, matrix, methodW) 
                                  for q in 1:size(matrix, 1) if q != o_o]), 
                       t_R, matrix, methodW)
    
    return (
        change_points = change_points,
        max_regret_alts = max_regret_alts,
        intervals = intervals,
        r_M = r_M,
        r_m = r_m
    )
end

# アベイルスペースに基づいて全ての代替案の順位変化点を分析
function analyze_all_alternatives_with_avail_space(utility::Matrix{Float64}, methodW, t_range::Tuple{Float64, Float64})
    n = size(utility, 1)
    matrix = create_minimax_R_Matrix(utility)
    
    # 各代替案の順位変化点を格納
    all_results = Dict{Int, Any}()
    
    for o_o in 1:n
        result = find_rank_change_points_from_tR(o_o, t_range[1], t_range[2], matrix, methodW)
        all_results[o_o] = result
    end
    
    # 全ての変化点をソート
    all_change_points = Float64[]
    for o_o in 1:n
        append!(all_change_points, all_results[o_o].change_points)
    end
    unique!(sort!(all_change_points))
    
    # アベイルスペースに基づく区間も収集
    all_intervals = Tuple{Float64, Float64}[]
    for o_o in 1:n
        append!(all_intervals, all_results[o_o].intervals)
    end
    unique!(all_intervals)
    
    # すべての重要な点（変化点と区間境界）をまとめる
    all_points = Float64[]
    append!(all_points, all_change_points)
    for (lower, upper) in all_intervals
        push!(all_points, lower)
        push!(all_points, upper)
    end
    unique!(sort!(all_points))
    
    # 各ポイントでの順位を計算
    rankings = Dict{Float64, Vector{Int}}()
    
    for t in all_points
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
        intervals = all_intervals,
        all_points = all_points,
        rankings = rankings
    )
end

end # module