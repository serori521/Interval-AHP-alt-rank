using Pkg
using DelimitedFiles,StatsBase
using CSV, DataFrames
using StringEncodings,Logging

#効用値行列を5行6列ごとに2次元配列として読み込み、配列に格納する関数
@inline function read_utility_value()
    # CSVファイルのパス
    csv_path = "/workspaces/Interval-AHP-alt-rank/data/効用値行列/u1/N=6_M=5/u.csv"
    # CSVファイルを読み込む（ヘッダーがないと指定）
    data = readdlm(csv_path, ',', Float64)
    df = DataFrame(data, :auto)
    # 5行6列ごとに2次元配列としてデータを格納
    utility_for_alta = []
    for i in 1:5:size(df, 1)
        push!(utility_for_alta, Matrix(df[i:i+4, :]))
    end
    return utility_for_alta
end

#真の区間重要度を読み込む関数
function read_true_weights(filename::String)
    csv_path = "/workspaces/Interval-AHP-alt-rank/data/Simp/a3/"*filename* "/Given_interval_weight.csv"
    data = readdlm(csv_path, ',', Float64)
    n = length(data)

    return (
        L = [data[i] for i in 1:2:n-1],
        R = [data[i] for i in 2:2:n]
    )
end

#それぞれの手法による区間重要度を読み込む関数
#repeat_num:1-1000 ,criteria_num:6 評価基準数 4-8
function read_method_weights(filename::String,repeat_num::Int,criteria_num::Int)
    csv_path = "/workspaces/Interval-AHP-alt-rank/data/Simp/a3/"*filename*"/Simp.csv"
    with_logger(NullLogger()) do
        data = readdlm(open(csv_path, enc"SHIFT_JIS"), ',', Float64; skipstart=3)
        # 各行ごとにタプルを作成
        result = [(
            L = data[i, 2:2:2+criteria_num*2-1],
            R = data[i, 3:2:2+criteria_num*2],
            adjacent = data[i, 2+criteria_num*2]
        ) for i in 1:repeat_num]

        return result

    end
end
