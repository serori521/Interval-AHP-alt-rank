include("calc_IPW.jl")

# 順位変化の結果をCSVに保存する関数
function save_rank_changes_to_csv(results, filename::String)
    # 結果を整形
    change_points = results.change_points
    rankings = results.rankings
    
    # すべての時点（初期点、変化点、終点）を取得
    t_points = sort(collect(keys(rankings)))
    
    # 各時点での順位を記録
    n = length(rankings[t_points[1]])  # 代替案の数
    
    # 出力用のデータフレームを作成
    output_data = DataFrame(
        t_point = Float64[],
        alt_1 = Int[],
        alt_2 = Int[],
        alt_3 = Int[],
        alt_4 = Int[],
        alt_5 = Int[]
    )
    
    for t in t_points
        # 各時点での順位を記録
        row = Dict("t_point" => t)
        for i in 1:n
            row["alt_$i"] = rankings[t][i]
        end
        push!(output_data, row)
    end
    
    # CSVに保存
    CSV.write(filename, output_data)
    println("Rank changes saved to $filename")
    
    return output_data
end
