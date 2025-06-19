module DecisionAnalysis

# export all key functions and types the user might need at the top level.

include("types.jl")
using .Types
export UtilityMatrix, MethodWeights, PairwiseRegretInput, SingleRegretResult,
    OverallRankingAtT, AnalysisResults, RankChangeDataPoint #, AlternativeSpecificAnalysis

include("io.jl")
using .DataIO
export read_utility_value, read_method_weights, read_true_weights, save_analysis_results_to_csv

include("utils.jl")
using .Utils
export find_optimal_trange

include("regret_calculations.jl")
using .RegretCalculations
export precompute_pairwise_inputs, calculate_single_pairwise_regret, get_max_regret_for_alternative

include("rank_analysis.jl")
using .RankAnalysis
export analyze_overall_rank_changes, get_overall_ranking_at_t # Exporting the main analysis function

end # module DecisionAnalysis