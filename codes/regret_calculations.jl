module RegretCalculations

using ..Types # Use .. to go up one level if files are in src/
using LinearAlgebra # For norm, dot, etc. if needed later, though not directly in this part

export precompute_pairwise_inputs, calculate_single_pairwise_regret, get_max_regret_for_alternative

# Precompute (U_q - U_o) and sorted indices for all pairs (o,q)
# Stores result for (o,q) where o is row, q is col.
function precompute_pairwise_inputs(utility_matrix::UtilityMatrix)::Matrix{Union{Nothing, PairwiseRegretInput}}
    n_alternatives = size(utility_matrix, 1)
    # cache[o, q] will store data for regret(o,q), i.e., based on (U_q - U_o)
    pairwise_inputs_cache = Matrix{Union{Nothing, PairwiseRegretInput}}(nothing, n_alternatives, n_alternatives)

    for o_idx in 1:n_alternatives
        for q_idx in 1:n_alternatives
            if o_idx == q_idx
                continue
            end
            # For regret of o being chosen over q (R(o,q)), we consider U_q - U_o
            diff_U_q_vs_o = utility_matrix[q_idx, :] - utility_matrix[o_idx, :]
            sorted_indices = sortperm(diff_U_q_vs_o, rev=true) # Descending order
            pairwise_inputs_cache[o_idx, q_idx] = PairwiseRegretInput(diff_U_q_vs_o, sorted_indices)
        end
    end
    return pairwise_inputs_cache
end

# Calculate regret R(o,q,t) = max { sum_j (u_qj - u_oj) y_j } 
# subject to y_j in [L_j*t, R_j*t] and sum(y_j) = t.
# diff_input = precomputed_pairwise_inputs[o_idx, q_idx]
function calculate_single_pairwise_regret(
    t::Float64,
    method_weights::MethodWeights,
    diff_input::PairwiseRegretInput
)::SingleRegretResult
    
    diff_U = diff_input.diff_U_q_vs_o
    sorted_diff_indices = diff_input.sorted_diff_indices

    Y_L_effective = method_weights.L .* t
    Y_R_effective = method_weights.R .* t
    width_Y_effective = Y_R_effective .- Y_L_effective
    num_criteria = length(method_weights.L)

    # Initial objective value with all y_j = L_j*t
    current_max_objective = sum(diff_U[j] * Y_L_effective[j] for j in 1:num_criteria)
    
    # Budget for increasing y_j beyond L_j*t, such that sum(y_j) = t
    # z_budget = t - sum(y_j_initial) = t - sum(L_j*t)
    z_budget = t - sum(Y_L_effective)

    # If sum(L_j*t) > t (i.e., sum(L_j) > 1), it's impossible to satisfy constraints.
    # This implies the max objective (regret) would be -Inf if truly maximizing under infeasibility.
    # Or, this 't' is outside the valid range where sum(y_j)=t is possible with y_j >= L_j*t.
    if z_budget < -1e-9 # Allow for small numerical errors
        # println("Warning: Negative z_budget ($z_budget) for t=$t. sum(L*t)=$(sum(Y_L_effective)). Infeasible.")
        return SingleRegretResult(-Inf, 0, 0.0)
    end
    
    intermediate_crit_idx_final = 0
    avail_space_final = 0.0

    # Iterate through criteria sorted by decreasing (u_qj - u_oj)
    # Allocate budget 'z_budget' to increase y_j from L_j*t up to R_j*t
    for k_loop in 1:num_criteria
        actual_criterion_idx = sorted_diff_indices[k_loop] # Index of the criterion in original L, R, diff_U
        
        # Max possible increase for this y_j is width_Y_effective[actual_criterion_idx]
        increase_amount = min(z_budget, width_Y_effective[actual_criterion_idx])
        
        current_max_objective += diff_U[actual_criterion_idx] * increase_amount
        z_budget -= increase_amount
        
        if abs(z_budget) < 1e-9 # Budget exhausted
            intermediate_crit_idx_final = actual_criterion_idx
            # Avail space is how much more y_j could have increased if budget was not limiting
            avail_space_final = width_Y_effective[actual_criterion_idx] - increase_amount
            break
        end
        # If loop finishes and z_budget is still > 0, it means sum(R_j*t) < t (sum(R_j) < 1)
        # This is also an infeasibility for sum(y_j)=t if y_j <= R_j*t. Max objective would be +Inf.
        if k_loop == num_criteria && z_budget > 1e-9
            # println("Warning: Positive z_budget ($z_budget) remaining at t=$t. sum(R*t)=$(sum(Y_R_effective)). Infeasible.")
            return SingleRegretResult(Inf, actual_criterion_idx, Inf) # Interm idx is the last one
        elseif k_loop == num_criteria # Budget perfectly used or slightly negative due to precision
            intermediate_crit_idx_final = actual_criterion_idx
            avail_space_final = 0.0 # No more budget, no more space in this criterion
        end
    end
    
    return SingleRegretResult(current_max_objective, intermediate_crit_idx_final, avail_space_final)
end


# For a given alternative 'o_idx', find its maximum regret against all other alternatives 'q_idx' at time 't'.
function get_max_regret_for_alternative(
    o_idx::Int,
    t::Float64,
    n_alternatives::Int,
    utility_matrix::UtilityMatrix, # Passed to calculate_single_pairwise_regret if not using full cache
    method_weights::MethodWeights,
    pairwise_inputs_cache::Matrix{Union{Nothing, PairwiseRegretInput}}
)::Tuple{Float64, Int} # Returns (max_regret_value, opponent_q_idx)

    max_r_val = -Inf
    opponent_idx_for_max_r = 0

    for q_idx in 1:n_alternatives
        if o_idx == q_idx
            continue
        end
        
        # R(o,q,t)
        regret_result = calculate_single_pairwise_regret(
            t,
            method_weights,
            pairwise_inputs_cache[o_idx, q_idx] # Pass the precomputed U_q - U_o
        )
        
        if regret_result.value > max_r_val
            max_r_val = regret_result.value
            opponent_idx_for_max_r = q_idx
        end
    end
    # If all regrets were -Inf (e.g., due to infeasibility), return that.
    return max_r_val, opponent_idx_for_max_r
end

end # module RegretCalculations