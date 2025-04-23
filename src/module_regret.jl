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
# 区間内での順位変化点とリグレット値の詳細な検出
function find_rank_change_in_interval(o_o::Int, t_lower::Float64, t_upper::Float64, matrix, methodW)
    n = size(matrix, 1)
    candidates = [q for q in 1:n if q != o_o]
    change_points = Float64[]
    max_regret_alts = Int[]
    
    # 現在の区間の左端を初期化
    current_t0 = t_lower
    
    # 左端での最大リグレットを与える代替案を特定
    regrets_at_t0 = [regret_value(o_o, q, current_t0, matrix, methodW) for q in candidates]
    p_idx = argmax(regrets_at_t0)
    o_p = candidates[p_idx]
    
    while current_t0 < t_upper
        # 右端でo_pよりも大きいリグレットを与える候補を特定
        regrets_at_t1 = [regret_value(o_o, q, t_upper, matrix, methodW) for q in candidates]
        O = [candidates[i] for i in 1:length(candidates) 
             if regrets_at_t1[i] > regret_value(o_o, o_p, t_upper, matrix, methodW)]
        
        # 候補がなければこの区間での変化点はなし（終了）
        if isempty(O)
            break
        end
        
        # 各候補との交点を計算
        crossing_points = Dict{Int, Float64}()
        for o_q in O
            t_cross = compute_crossing_point(o_o, o_p, o_q, current_t0, t_upper, matrix, methodW)
            if t_cross !== nothing && t_cross > current_t0 && t_cross <= t_upper
                crossing_points[o_q] = t_cross
            end
        end
        
        # 交点がなければ終了
        if isempty(crossing_points)
            break
        end
        
        # 最小交点とその代替案を特定
        o_r, t_min = findmin(crossing_points)
        
        # 変化点を記録
        push!(change_points, t_min)
        push!(max_regret_alts, o_r)
        
        # 更新して繰り返し
        current_t0 = t_min
        o_p = o_r
    end
    
    return change_points, max_regret_alts
end

# アベイルスペースに基づいてt^Rから順に区間を分析
function find_rank_change_points_from_tR(o_o::Int, t_L::Float64, t_R::Float64, matrix, methodW)
    # 初期化
    all_change_points = Float64[]
    all_max_regret_alts = Int[]
    all_max_regret_values = Float64[]
    intervals = Tuple{Float64, Float64}[]
    interval_data = Vector{Any}()
    
    # 現在のt値を初期化
    current_t = t_R
    
    while current_t > t_L
        # 次のt値を計算
        next_t, min_availspace = compute_next_t(current_t, matrix, methodW)
        
        # 下限t_Lを下回らないようにする
        next_t = max(next_t, t_L)
        
        # 区間を記録
        push!(intervals, (next_t, current_t))
        
        # 区間内の順位変化点を検出
        interval_result = find_rank_change_in_interval(o_o, current_t, next_t, matrix, methodW)
        push!(interval_data, interval_result)
        println(interval_result)
        # 全体の結果に追加
        append!(all_change_points, interval_result.change_points)
        append!(all_max_regret_alts, interval_result.max_regret_alts)
        append!(all_max_regret_values, interval_result.max_regret_values)
        
        # 次の区間へ
        current_t = next_t
        
        # t_Lに到達したら終了
        if isapprox(current_t, t_L, atol=1e-8)
            break
        end
    end
    
    # 全体的な最大・最小リグレット値
    r_M = maximum([interval.r_max for interval in interval_data])
    r_m = minimum([interval.r_min for interval in interval_data])
    
    return (
        change_points = all_change_points,
        max_regret_alts = all_max_regret_alts,
        max_regret_values = all_max_regret_values,
        intervals = intervals,
        interval_data = interval_data,
        r_M = r_M,
        r_m = r_m
    )
end

# 代替案間の関係を分析し、全体的な順位変化を検出
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
    
    # 代替案間の比較による追加の順位変化点を検出
    additional_change_points = Float64[]
    
    # 各代替案ペアの区間でリグレット範囲が重なる部分を検出
    for i in 1:n
        for j in i+1:n
            # 各代替案の区間データを取得
            i_intervals = all_results[i].interval_data
            j_intervals = all_results[j].interval_data
            
            # 各区間の組み合わせでリグレット範囲の重なりを確認
            for i_data in i_intervals
                for j_data in j_intervals
                    # リグレット範囲が重なる場合
                    if !(i_data.r_max < j_data.r_min || i_data.r_min > j_data.r_max)
                        # 線形関数のパラメータを計算し、交点を求める
                        # (簡略化のため、詳細な実装は省略)
                        # この部分は代替案間の比較による順位変化点の検出となります
                    end
                end
            end
        end
    end
    
    # 追加の変化点を統合
    append!(all_change_points, additional_change_points)
    unique!(sort!(all_change_points))
    if change_points == []
        all_change_points = [0]
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