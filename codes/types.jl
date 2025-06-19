module Types

export UtilityMatrix, MethodWeights, PairwiseRegretInput, SingleRegretResult, OverallRankingAtT, RankChangeDataPoint, AlternativeSpecificAnalysis, AnalysisResults

# Matrix of utility values (alternatives x criteria)
const UtilityMatrix = Matrix{Float64}

# Interval weights for a method
# type.jl内で定義
const MethodWeights = @NamedTuple{L::Vector{Float64}, R::Vector{Float64}, adjacent::Float64}

# Precomputed difference U_q - U_o and sorted indices for a pair (o,q)
struct PairwiseRegretInput
    diff_U_q_vs_o::Vector{Float64} # utility_qj - utility_oj for all criteria j
    sorted_diff_indices::Vector{Int}    # Indices of criteria, sorted by diff_U_q_vs_o descending
end

# Result of calculating regret for a single pair (o,q) at a specific t
struct SingleRegretResult
    value::Float64
    intermediate_criterion_idx::Int # The criterion index that exhausted the budget 'z'
    avail_space::Float64            # Remaining capacity in the intermediate_criterion's weight interval
end

# Overall ranking at a specific t value
struct OverallRankingAtT
    t::Float64
    ranking::Vector{Int}              # Alternative indices sorted by min max_regret
    max_regrets::Vector{Float64}      # Max regret for each alternative at this t
    # contributing_opponents::Vector{Int} # For each alt, which opponent gave the max regret
end

# Structure to hold all rank change analysis results
struct RankChangeDataPoint
    t_value::Float64
    ranking::Vector{Int}
    # max_regrets::Vector{Float64} # Optional: store max regrets at these specific change points
end

struct AlternativeSpecificAnalysis
    change_points::Vector{Float64} # t-values where max regret opponent for this alt changes
    # ... other specific data if needed ...
end

struct AnalysisResults
    utility_matrix::UtilityMatrix
    method_weights::MethodWeights
    t_range_analyzed::Tuple{Float64,Float64}

    # Points where overall ranking actually changes due to intersections
    true_overall_rank_change_points::Vector{RankChangeDataPoint}

    # Rankings at all evaluated critical t-points (for plotting segments)
    all_evaluated_rankings::Vector{OverallRankingAtT}

    # Detailed analysis for each alternative (optional, if still needed for diagnostics)
    # results_per_alternative::Dict{Int, AlternativeSpecificAnalysis}

    # Intervals derived from avail_space logic (used to define linearity of regret functions)
    # avail_space_intervals::Vector{Tuple{Float64, Float64}}
end

end # module Types