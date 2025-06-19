module DataIO

using DelimitedFiles, DataFrames, StringEncodings
# Add CSV if you prefer it for writing DataFrames, otherwise DelimitedFiles can handle it.
# import Pkg; Pkg.add("CSV"); using CSV 
using ..Types # Assuming Types.jl is in the same directory level

export read_utility_value, read_method_weights, read_true_weights, save_analysis_results_to_csv

function read_utility_value(filepath::String="/workspaces/inulab_julia_devcontainer/data/効用値行列/u1/N=6_M=5/u.csv")::Vector{Matrix{Float64}}
    data = readdlm(filepath, ',', Float64)
    # Assuming 5 alternatives (rows per matrix)
    num_alternatives = 5
    utility_matrices = UtilityMatrix[]

    num_criteria = size(data, 2)

    for i in 1:num_alternatives:size(data, 1)
        push!(utility_matrices, data[i:i+num_alternatives-1, 1:num_criteria])
    end
    return utility_matrices # Or utility_matrices[1] if only one is used at a time
end

function read_method_weights(filename::String, repeat_num::Int, criteria_num::Int)
    csv_path = "/workspaces/inulab_julia_devcontainer/data/Simp/a3/" * filename * "/Simp.csv"

    # StringEncodingsを直接使用
    io = open(csv_path, enc"SHIFT_JIS", "r")
    data = readdlm(io, ',', Float64; skipstart=3)
    close(io)

    # メモリ効率の良いデータ生成
    result = Vector{NamedTuple}(undef, repeat_num)
    for i in 1:repeat_num
        result[i] = (
            L=data[i, 2:2:2+criteria_num*2-1],
            R=data[i, 3:2:2+criteria_num*2],
            adjacent=data[i, 2+criteria_num*2]
        )
    end
    return result
end

function read_true_weights(filepath::String)::MethodWeights
    data = readdlm(filepath, ',', Float64)
    n = length(data)
    L = [data[i] for i in 1:2:n-1]
    R = [data[i] for i in 2:2:n]
    return MethodWeights(L, R) #, nothing) # Assuming 'adjacent' is not in this file
end

function save_analysis_results_to_csv(results::AnalysisResults, filename::String="rank_change_details.csv")
    rank_data = results.true_overall_rank_change_points
    if isempty(rank_data)
        println("No rank change data to save for $filename.")
        return DataFrame()
    end

    n_alternatives = length(rank_data[1].ranking)

    # Prepare data for DataFrame
    t_points_col = [rd.t_value for rd in rank_data]

    df = DataFrame(t_change_point=t_points_col)
    for alt_id in 1:n_alternatives
        alt_rank_col = Symbol("alt_$(alt_id)_rank")
        df[!, alt_rank_col] = Int[] # Initialize column
    end

    for rd in rank_data
        # rd.ranking is like [alt_at_rank1, alt_at_rank2, ...]
        # We need to find the rank of alt_1, rank of alt_2, etc.
        current_alt_ranks = zeros(Int, n_alternatives)
        for rank_pos in 1:length(rd.ranking)
            alt_in_this_pos = rd.ranking[rank_pos]
            if 1 <= alt_in_this_pos <= n_alternatives
                current_alt_ranks[alt_in_this_pos] = rank_pos
            end
        end

        new_row = (t_change_point=rd.t_value,)
        for alt_id in 1:n_alternatives
            new_row = merge(new_row, (Symbol("alt_$(alt_id)_rank") => current_alt_ranks[alt_id],))
        end
        push!(df, new_row, promote=true) # promote needed if initial df was empty with specific types
    end

    # Ensure t_change_point is the first column if it got reordered
    cols_ordered = [:t_change_point]
    for alt_id in 1:n_alternatives
        push!(cols_ordered, Symbol("alt_$(alt_id)_rank"))
    end
    df_ordered = df[!, cols_ordered]

    try

        header_str = join(names(df_ordered), ',')
        open(filename, "w") do io
            println(io, header_str)
            writedlm(io, Iterators.eachrow(df_ordered), ',')
        end
        println("Analysis results saved to $filename")
    catch e
        println("Error saving CSV to $filename: $e")
    end
    return df_ordered
end

end # module DataIO