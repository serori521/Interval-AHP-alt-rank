using Pkg
using DelimitedFiles,StatsBase
using CSV, DataFrames
using StringEncodings

#効用値行列を5行6列ごとに2次元配列として読み込み、配列に格納する関数
@inline function read_utility_value()
    # CSVファイルのパス
    csv_path = "/workspaces/Interval-AHP-alt-rank/data/効用値行列/u1/N=6_M=5/u.csv"
    # CSVファイルを読み込む（ヘッダーがないと指定）
    df = CSV.File(csv_path, header=false) |> DataFrame
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
    data = first(CSV.File(csv_path, header=false))
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
        data = CSV.File(open(csv_path, enc"SHIFT_JIS"); header=1:3, delim=',') |> DataFrame
        # B列の4行目を起点に、指定された列数と行数のデータを抽出
        values = data[1:repeat_num, 2:1+criteria_num*2]
        #　adjacent は各区間の幅の総和を表す
        adjacent_values = data[1:repeat_num, 2+criteria_num*2]
        # 交互にLとRに格納（行列形式repeat行criteria列）
        L = reshape([values[row, col] for col in 1:2:criteria_num*2-1 for row in 1:repeat_num], repeat_num, criteria_num)
        R = reshape([values[row, col] for col in 2:2:criteria_num*2 for row in 1:repeat_num], repeat_num, criteria_num)
        #返り値はタプル
        return (
            L = L,
            R = R,
            adjacent = Vector(adjacent_values)
        );
    end
end