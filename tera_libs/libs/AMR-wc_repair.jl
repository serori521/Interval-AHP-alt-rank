"""
提案法eMMRw/c, gMMRw/c, lMMRw/c
MMR_kai(PCM, method)とすることで，区間重要度が求められる．
"""

using IntervalArithmetic
using JuMP
import HiGHS

using Plots
include("./crisp-pcm.jl")
include("./nearly-equal.jl")
include("./solve-deterministic-ahp.jl")


MMRE_Individual_kai = @NamedTuple{
    # 区間重みベクトル
    s::T,
    centers::Matrix{T},
    l::Matrix{T},
    wᴸ::Vector{T}, wᵁ::Vector{T},
    W::Vector{Interval{T}} # ([wᵢᴸ, wᵢᵁ])
} where {T <: Real}

# 任意の行と列を削除
@inline function remove_row_col(A::Matrix{T}, row::Int, col::Int)::Matrix{T} where {T <: Real}
    m, n = size(A)

    # 行を除外
    new_matrix = A[setdiff(1:m, row), :]
    # 列を除外
    result_matrix = new_matrix[:, setdiff(1:n, col)]
        
    return result_matrix
end

@inline function phase1_kai(A::Matrix{T}, method::Function)::Matrix{T} where {T <: Real}

    m, n = size(A)

    # Phase 1
    Wᶜ = Matrix{T}(undef, m, n) 
    for k = 1:n
        removed_matrix = remove_row_col(A, k, k)
        W = method(removed_matrix)

        # Wᶜₖを1として挿入
        insert!(W, k, 1.0)
        Wᶜ[:, k] = W
    end

    return Wᶜ
end

# Phase2のループの中の部分
@inline function phase2_jump_kai(A::Matrix{T}, Wᶜ::Matrix{T}, k::Int, n::Int)::T where {T <: Real}
    ε = 1e-6 # << 1

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, L[i=1:n] ≥ 0) # ここはイプシロンを０に変更
        @variable(model, t ≥ 1+ε) # 極限をとるのでOK
        Lₖ = L[k]
        
        for j = filter(j -> j != k, 1:n)
            for i = filter(i -> i != j && i != k, 1:n)
                aᵢⱼ = A[i,j]
                wᵢᶜ = Wᶜ[i,k]; wⱼᶜ = Wᶜ[j,k]
                Lᵢ = L[i]; Lⱼ = L[j]
                @constraint(model, aᵢⱼ*(wⱼᶜ-Lⱼ) ≤ wᵢᶜ+Lᵢ)
            end
        end

        for j = filter(j -> j != k, 1:n)
            aₖⱼ = A[k, j]
            Lⱼ = L[j]
            wⱼᶜ = Wᶜ[j,k]
            @constraint(model, aₖⱼ*(wⱼᶜ-Lⱼ) ≤ t-1 +Lₖ)
        end

        for i = filter(i -> i != k, 1:n)
            aᵢₖ = A[i,k]
            Lᵢ = L[i]
            wᵢᶜ = Wᶜ[i,k]
            @constraint(model, aᵢₖ*(t-1-Lₖ) ≤ wᵢᶜ+Lᵢ)
            @constraint(model, wᵢᶜ - Lᵢ ≥ t*ε)
        end

        for j = filter(j -> j != k, 1:n)
            Lⱼ = L[j];
            wⱼᶜ = Wᶜ[j,k];
            # Sᵁ = Σ(μₖ*wᵢᶜ+lᵢ)
            # Sᴸ = Σ(μₖ*wᵢᶜ-lᵢ)
            Sᵁ = sum(map(i -> Wᶜ[i,k]+L[i], filter(i -> i != j && i != k, 1:n)))
            Sᴸ = sum(map(i -> Wᶜ[i,k]-L[i], filter(i -> i != j && i != k, 1:n)))
            @constraint(model, (t-1)+Lₖ + Sᵁ + wⱼᶜ-Lⱼ ≥ t)
            @constraint(model, (t-1)-Lₖ + Sᴸ + wⱼᶜ+Lⱼ ≤ t)
        end

        # Sᵁ_dash = Σ(μₖ*wᵢᶜ+lᵢ)
        # Sᴸ_dash = Σ(μₖ*wᵢᶜ-lᵢ)
        Sᵁ_dash = sum(map(i -> Wᶜ[i,k]+L[i], filter(i -> i != k, 1:n)))
        Sᴸ_dash = sum(map(i -> Wᶜ[i,k]-L[i], filter(i -> i != k, 1:n)))
        @constraint(model, Sᵁ_dash + (t-1)-Lₖ ≥ t)
        @constraint(model, Sᴸ_dash + (t-1)+Lₖ ≤ t)
        @constraint(model, (t-1)-Lₖ ≥ t*ε)

        dₖ  = sum(map(j -> L[j], filter(j -> j != k, 1:n)))
        @objective(model, Min, dₖ)

        optimize!(model)

        dₖ⃰ = sum(map(j -> value.(L[j]), filter(j -> j != k, 1:n)))

        # for i = 1:n
        #     println("k=$k, l[$i] = ", value.(l[i]))
        # end
        # println("k=$k, μₖ = ", value.(μₖ))
        # println("k=$k, dₖ⃰ = ", dₖ⃰)
        
        return dₖ⃰
        
    finally
        empty!(model)
    end
end

# Phase3の戻り値
# 変更後
phase3_jump_result_kai = @NamedTuple{
    μₖ⃰::T, l⃰::Vector{T},
} where {T <: Real}

@inline function phase3_jump_kai(A::Matrix{T}, Wᶜ::Matrix{T}, d⃰::T, k::Int, n::Int)::phase3_jump_result_kai{T} where {T <: Real}
    ε = 1e-6
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, l[i=1:n] ≥ ε)
        @variable(model, ε ≤ μₖ ≤ 1-ε)
        lₖ = l[k]
        
        # 制約1: aᵢⱼ(μₖwⱼᶜ - lⱼ) ≤ μₖwᵢᶜ + lᵢ
        for j = filter(j -> j != k, 1:n)
            for i = filter(i -> i != j && i != k, 1:n)
                aᵢⱼ = A[i,j]
                wᵢᶜ = Wᶜ[i,k]; wⱼᶜ = Wᶜ[j,k]
                lᵢ = l[i]; lⱼ = l[j]
                @constraint(model, aᵢⱼ*(μₖ*wⱼᶜ-lⱼ) ≤ μₖ*wᵢᶜ+lᵢ)
            end
        end

        # 制約2: aₖⱼ(μₖwⱼᶜ - lⱼ) ≤ 1-μₖ + lₖ
        for j = filter(j -> j != k, 1:n)
            aₖⱼ = A[k, j]
            lⱼ = l[j]
            wⱼᶜ = Wᶜ[j,k]
            @constraint(model, aₖⱼ*(μₖ*wⱼᶜ-lⱼ) ≤ 1-μₖ +lₖ)
        end

        # 制約3: aᵢₖ((1-μₖ) - lₖ) ≤ μₖwᵢᶜ + lᵢ
        for i = filter(i -> i != k, 1:n)
            aᵢₖ = A[i,k]
            lᵢ = l[i]
            wᵢᶜ = Wᶜ[i,k]
            @constraint(model, aᵢₖ*(1-μₖ-lₖ) ≤ μₖ*wᵢᶜ+lᵢ)
            @constraint(model, μₖ*wᵢᶜ - lᵢ ≥ ε)
        end

        # 制約4,5: 和の制約
        for j = filter(j -> j != k, 1:n)
            lⱼ = l[j]
            wⱼᶜ = Wᶜ[j,k]
            Sᵁ = sum(map(i -> μₖ*Wᶜ[i,k]+l[i], filter(i -> i != j && i != k, 1:n)))
            Sᴸ = sum(map(i -> μₖ*Wᶜ[i,k]-l[i], filter(i -> i != j && i != k, 1:n)))
            @constraint(model, (1-μₖ)+lₖ + Sᵁ + μₖ*wⱼᶜ-lⱼ ≥ 1)
            @constraint(model, (1-μₖ)-lₖ + Sᴸ + μₖ*wⱼᶜ+lⱼ ≤ 1)
        end

        # 制約6,7: 全体の和の制約
        Sᵁ_dash = sum(map(i -> μₖ*Wᶜ[i,k]+l[i], filter(i -> i != k, 1:n)))
        Sᴸ_dash = sum(map(i -> μₖ*Wᶜ[i,k]-l[i], filter(i -> i != k, 1:n)))
        @constraint(model, Sᵁ_dash + (1-μₖ)-lₖ ≥ 1)
        @constraint(model, Sᴸ_dash + (1-μₖ)+lₖ ≤ 1)
        
        # 制約8: Σlⱼ = d⃰
        Σl = sum(map(j-> l[j], filter(j -> j != k, 1:n)))
        @constraint(model, Σl == d⃰)  # 等式制約に変更
        
        # 制約9: 非負制約
        @constraint(model, (1-μₖ)-lₖ ≥ ε)

        @objective(model, Min, sum(l))

        optimize!(model)
        μₖ⃰ = value(μₖ)
        l⃰ = value.(l)
        return (μₖ⃰ = μₖ⃰, l⃰ = l⃰)

    finally
        empty!(model)
    end
end

# 提案手法 MMR-E, MMR-G, MMR-A
@inline function MMR_kai(A::Matrix{T}, method::Function)::MMRE_Individual_kai{T} where {T <: Real}

    if !isCrispPCM(A)
        throw(ArgumentError("A is not a crisp PCM"))
    end

    m, n = size(A)

    # Phase 1
    Wᶜ = phase1_kai(A, method)

    # Phase 2
    d⃰ = Vector{T}(undef, n) 
    for k = 1:n
        d⃰[k] = phase2_jump_kai(A, Wᶜ, k, n)
    end

    # MMR_kai関数内のPhase 3部分
    for k = 1:n
        (μₖ⃰, l⃰) = phase3_jump_kai(A, Wᶜ, d⃰[k], k, n)  # 戻り値の受け取り方を変更

        for i = 1:n
            lᵢ⃰ = l⃰[i]
            wᵢᶜ = Wᶜ[i,k]
            if i==k
                wᴸᵢ = (1-μₖ⃰) - lᵢ⃰
                wᵁᵢ = (1-μₖ⃰) + lᵢ⃰
                centers[i, k] = 1-μₖ⃰
            else
                wᴸᵢ = μₖ⃰*wᵢᶜ - lᵢ⃰
                wᵁᵢ = μₖ⃰*wᵢᶜ + lᵢ⃰
                centers[i, k] = μₖ⃰*wᵢᶜ
            end
            l[i, k] = lᵢ⃰
            wᴸ[i, k] = wᴸᵢ
            wᵁ[i, k] = wᵁᵢ
        end
    end


    # Phase 4
    w̅̅ᴸ = Vector{T}(undef, n)
    w̅̅ᵁ = Vector{T}(undef, n)
    W̅̅ = Vector{Interval{T}}(undef, n)

    for i = 1:n
        w̅̅ᴸ[i] = mean(wᴸ[i, :])
        w̅̅ᵁ[i] = mean(wᵁ[i, :])

        # precision error 対応
        if w̅̅ᴸ[i] > w̅̅ᵁ[i]
            w̅̅ᴸ[i] = w̅̅ᵁ[i]
        end

    end

    for i = 1:n
        W̅̅[i] = interval(w̅̅ᴸ[i], w̅̅ᵁ[i])
    end

    return (
        wᴸ=w̅̅ᴸ, wᵁ=w̅̅ᵁ,
        W=W̅̅
    )

end