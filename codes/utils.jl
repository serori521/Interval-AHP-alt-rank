module Utils

using ..Types

export find_optimal_trange

# From calc_IPW.jl
function find_optimal_trange(L::Vector{Float64}, R::Vector{Float64})
    n = length(L)

    max_denominator_val = -Inf
    for j in 1:n
        current_sum = sum(L[i] for i in 1:n if i != j) + R[j]
        max_denominator_val = max(max_denominator_val, current_sum)
    end
    t_range_start = 1 / max_denominator_val # This was t_R in original, maps to t_range[1]
    min_denominator_val = Inf
    for j in 1:n
        current_sum = sum(R[i] for i in 1:n if i != j) + L[j]
        min_denominator_val = min(min_denominator_val, current_sum)
    end
    t_range_end = 1 / min_denominator_val # This was t_L in original, maps to t_range[2]

    # Ensure t_range_start < t_range_end
    return min(t_range_start, t_range_end), max(t_range_start, t_range_end)
end

end # module Utils