module RankAnalysis

using ..Types, ..RegretCalculations, ..Utils
# For root finding, consider adding Roots.jl to your Project.toml
# import Pkg; Pkg.add("Roots"); using Roots 
using Printf # For debugging output

export analyze_overall_rank_changes

# Helper to get overall ranking at a specific t
function get_overall_ranking_at_t(
    t::Float64,
    n_alternatives::Int,
    utility_matrix::UtilityMatrix,
    method_weights::MethodWeights,
    pairwise_inputs_cache::Matrix{Union{Nothing,PairwiseRegretInput}}
)::OverallRankingAtT

    max_regrets_for_alts = Vector{Float64}(undef, n_alternatives)
    # contributing_opponents_for_alts = Vector{Int}(undef, n_alternatives) # If needed

    for o_idx in 1:n_alternatives
        max_r, _ = get_max_regret_for_alternative(
            o_idx, t, n_alternatives, utility_matrix, method_weights, pairwise_inputs_cache
        )
        max_regrets_for_alts[o_idx] = max_r
        # contributing_opponents_for_alts[o_idx] = opponent_q_idx
    end

    # Sort alternatives by their max_regret values (ascending)
    # sortperm gives indices that would sort the array
    ranking_indices = sortperm(max_regrets_for_alts)

    return OverallRankingAtT(t, ranking_indices, max_regrets_for_alts[ranking_indices])
end


# Function to find intervals where a single pairwise regret R(o,q,t) is linear
# This is based on changes to intermediate_criterion_idx or avail_space
function get_linearity_intervals_for_pairwise_regret(
    o_idx::Int, q_idx::Int,
    t_start::Float64, t_end::Float64,
    method_weights::MethodWeights,
    diff_input_oq::PairwiseRegretInput
)
    intervals = Tuple{Float64,Float64,Int,Float64}[] # (t_low, t_high, interm_idx, initial_slope_approximation)

    current_t = t_start
    tol = 1e-9 # Tolerance for t comparisons

    while current_t < t_end - tol
        res_current = calculate_single_pairwise_regret(current_t, method_weights, diff_input_oq)

        # Estimate where the current linearity might break based on avail_space
        # next_t_break_estimate needs more sophisticated calculation based on how avail_space affects linearity.
        # For now, we simplify: the avail_space logic primarily defines the piece-wise nature.
        # The critical points are where interm_idx changes.
        # A simpler approach for now: assume linearity between points where interm_idx is known to change.
        # This needs the global set of critical points first.

        # This function is hard to implement without a global list of points where interm_idx might change.
        # Let's defer its full implementation or simplify the main algorithm's dependency on it.

        # For now, let's assume that the main `analyze_overall_rank_changes` will provide small enough
        # intervals `(t_a, t_b)` within which `intermediate_criterion_idx` for relevant pairs is constant.
        # If not, this function would be needed to find those sub-intervals.
        # For simplicity, we'll return the whole interval, assuming linearity for now.
        # The root finding will then operate on this assumed linear segment.
        # A more robust solution would actually find these break points.

        # Simplified: we'll rely on the global critical points to define small intervals.
        # This function is more about getting A, B for R(o,q,t) = A*t + B in an interval.

        res_end_interval = calculate_single_pairwise_regret(t_end, method_weights, diff_input_oq)

        if abs(t_end - current_t) < tol # Interval too small
            push!(intervals, (current_t, t_end, res_current.intermediate_criterion_idx, 0.0)) # Slope 0 for tiny interval
            break
        end

        # Check if intermediate criterion is the same at start and end of the small interval
        # If so, assume linearity. This is an approximation if the interval is not small enough.
        if res_current.intermediate_criterion_idx == res_end_interval.intermediate_criterion_idx || res_current.value == res_end_interval.value
            slope = (res_end_interval.value - res_current.value) / (t_end - current_t)
            push!(intervals, (current_t, t_end, res_current.intermediate_criterion_idx, slope))
        else
            # Intermediate criterion changed, split the interval (e.g., midpoint) and recurse or refine
            # This requires a more complex iterative approach.
            # For now, we'll still approximate with a single line.
            # A better way: the main function must ensure intervals are small enough.
            slope = (res_end_interval.value - res_current.value) / (t_end - current_t)
            push!(intervals, (current_t, t_end, res_current.intermediate_criterion_idx, slope)) # Or mark as non-linear
        end
        current_t = t_end # Move to next interval (this simplified version only does one)
    end
    if isempty(intervals) && t_start ≈ t_end # Single point
        res_current = calculate_single_pairwise_regret(t_start, method_weights, diff_input_oq)
        push!(intervals, (t_start, t_end, res_current.intermediate_criterion_idx, 0.0))
    end

    return intervals
end


function analyze_overall_rank_changes(
    utility_matrix::UtilityMatrix,
    method_weights::@NamedTuple{L::Vector{Float64}, R::Vector{Float64}, adjacent::Float64},
    t_range::Tuple{Float64,Float64}
)::AnalysisResults

    n_alternatives = size(utility_matrix, 1)
    pairwise_inputs_cache = precompute_pairwise_inputs(utility_matrix)

    # --- Step 1: Collect all potential critical t-points ---
    critical_t_set = Set{Float64}()
    push!(critical_t_set, t_range[1], t_range[2])

    for i in 0:100 # e.g., 100 steps, can be adaptive
        push!(critical_t_set, t_range[1] + i / 100 * (t_range[2] - t_range[1]))
    end

    # Add avail_space boundaries (iterative refinement)
    # This simulates the original approach's interval generation
    temp_t_points = sort(collect(critical_t_set))
    new_critical_points_found = true
    iter_count = 0
    max_iters = 5 # Limit iterations to prevent excessive points

    while new_critical_points_found && iter_count < max_iters
        new_critical_points_found = false
        iter_count += 1
        current_t_idx = 1
        while current_t_idx <= length(temp_t_points)
            t_curr = temp_t_points[current_t_idx]
            min_delta_t_overall = Inf

            for o_idx in 1:n_alternatives
                for q_idx in 1:n_alternatives
                    if o_idx == q_idx
                        continue
                    end
                    regret_res = calculate_single_pairwise_regret(t_curr, method_weights, pairwise_inputs_cache[o_idx, q_idx])
                    if regret_res.value == -Inf || regret_res.value == Inf || regret_res.avail_space == Inf || regret_res.avail_space <= 1e-9
                        continue
                    end

                end
            end
            # For now, this part of critical point generation is simplified.
            # True intersections are the most important.
            current_t_idx += 1
        end
    end


    # Add intersections of R(o, q1, t) and R(o, q2, t) for each o
    for o_idx in 1:n_alternatives
        for q1_idx in 1:n_alternatives
            if q1_idx == o_idx
                continue
            end
            for q2_idx in (q1_idx+1):n_alternatives # Avoid redundant pairs and self-comparison
                if q2_idx == o_idx
                    continue
                end

                # Define h(t) = R(o,q1,t) - R(o,q2,t)
                h = t_val -> calculate_single_pairwise_regret(t_val, method_weights, pairwise_inputs_cache[o_idx, q1_idx]).value -
                             calculate_single_pairwise_regret(t_val, method_weights, pairwise_inputs_cache[o_idx, q2_idx]).value

            end
        end
    end

    sorted_critical_t_points = sort(collect(filter(p -> t_range[1] <= p <= t_range[2], critical_t_set)))
    unique!(sorted_critical_t_points) # Ensure uniqueness and sorted

    # --- Step 2: Evaluate rankings and find true change points ---
    all_evaluated_rankings_list = OverallRankingAtT[]
    true_rank_change_data_points = RankChangeDataPoint[]

    if isempty(sorted_critical_t_points)
        # Ensure at least start and end of t_range are evaluated if critical set is empty
        push!(sorted_critical_t_points, t_range[1])
        if t_range[1] != t_range[2]
            push!(sorted_critical_t_points, t_range[2])
        end
        unique!(sorted_critical_t_points)
    end

    # Initial ranking at the start of the range
    prev_ranking_data = get_overall_ranking_at_t(sorted_critical_t_points[1], n_alternatives, utility_matrix, method_weights, pairwise_inputs_cache)
    push!(all_evaluated_rankings_list, prev_ranking_data)
    push!(true_rank_change_data_points, RankChangeDataPoint(prev_ranking_data.t, prev_ranking_data.ranking))

    for i in 2:length(sorted_critical_t_points)
        t_prev = sorted_critical_t_points[i-1]
        t_curr = sorted_critical_t_points[i]

        # Ranking at t_curr (use previous if t_curr is too close to t_prev)
        local curr_ranking_data::OverallRankingAtT
        if abs(t_curr - t_prev) < 1e-9
            curr_ranking_data = prev_ranking_data
        else
            curr_ranking_data = get_overall_ranking_at_t(t_curr, n_alternatives, utility_matrix, method_weights, pairwise_inputs_cache)
        end
        push!(all_evaluated_rankings_list, curr_ranking_data)

        if curr_ranking_data.ranking != prev_ranking_data.ranking


            actual_change_t_in_interval = t_curr # Default to t_curr


            # If the last recorded change point is not t_curr (to avoid duplicates if multiple critical points yield same ranking then change)
            if isempty(true_rank_change_data_points) || last(true_rank_change_data_points).t_value != t_curr
                push!(true_rank_change_data_points, RankChangeDataPoint(t_curr, curr_ranking_data.ranking))
            end
        end
        prev_ranking_data = curr_ranking_data
    end

    return AnalysisResults(
        utility_matrix,
        method_weights,
        t_range,
        true_rank_change_data_points,
        all_evaluated_rankings_list
    )
end

end # module RankAnalysis