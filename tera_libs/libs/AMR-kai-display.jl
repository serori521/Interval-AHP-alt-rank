"""
提案法eAMRw/c, gAMRw/c, lAMRw/c
各モデルの最適化結果観察用
"""

using IntervalArithmetic
using JuMP
import HiGHS
using Statistics

using Plots
include("./crisp-pcm.jl")
include("./nearly-equal.jl")
include("./solve-deterministic-ahp.jl")


AMRkai_Individual = @NamedTuple{
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

@inline function AMRkai_phase1(A::Matrix{T}, method::Function)::Matrix{T} where {T <: Real}

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
@inline function AMRkai_phase2_jump(A::Matrix{T}, Wᶜ::Matrix{T}, k::Int, n::Int)::T where {T <: Real}
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
AMRkai_phase3_jump_result = @NamedTuple{
    # 区間重みベクトル
    tₖ⃰::T, L⃰::Vector{T},
} where {T <: Real}

# Phase3のループの中の部分
@inline function AMRkai_phase3_jump(A::Matrix{T}, Wᶜ::Matrix{T}, d⃰::T, k::Int, n::Int)::AMRkai_phase3_jump_result{T} where {T <: Real}
    ε = 1e-6 # << 1

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, L[i=1:n] ≥ 0)
        @variable(model, t ≥ 1+ε)
        Lₖ = L[k]
        
        for j = filter(j -> j != k, 1:n)
            for i = filter(i -> i != j && i != k, 1:n)
                aᵢⱼ = A[i,j]
                wᵢᶜ = Wᶜ[i,k]; wⱼᶜ = Wᶜ[j,k]
                Lᵢ = L[i]; Lⱼ = L[j]
                @constraint(model, aᵢⱼ*((t-1)*wⱼᶜ-Lⱼ) ≤ (t-1)*wᵢᶜ+Lᵢ)
            end
        end

        for j = filter(j -> j != k, 1:n)
            aₖⱼ = A[k, j]
            Lⱼ = L[j]
            wⱼᶜ = Wᶜ[j,k]
            @constraint(model, aₖⱼ*((t-1)*wⱼᶜ-Lⱼ) ≤ 1+Lₖ)
        end

        for i = filter(i -> i != k, 1:n)
            aᵢₖ = A[i,k]
            Lᵢ = L[i]
            wᵢᶜ = Wᶜ[i,k]
            @constraint(model, aᵢₖ*(1-Lₖ) ≤ (t-1)*wᵢᶜ+Lᵢ)
            @constraint(model, (t-1)*wᵢᶜ - Lᵢ ≥ t*ε)
        end

        for j = filter(j -> j != k, 1:n)
            Lⱼ = L[j];
            wⱼᶜ = Wᶜ[j,k];
            # Sᵁ = Σ(μₖ*wᵢᶜ+lᵢ)
            # Sᴸ = Σ(μₖ*wᵢᶜ-lᵢ)
            Sᵁ = sum(map(i -> (t-1)*Wᶜ[i,k]+L[i], filter(i -> i != j && i != k, 1:n)))
            Sᴸ = sum(map(i -> (t-1)*Wᶜ[i,k]-L[i], filter(i -> i != j && i != k, 1:n)))
            @constraint(model, 1+Lₖ + Sᵁ + (t-1)*wⱼᶜ-Lⱼ ≥ t)
            @constraint(model, 1-Lₖ + Sᴸ + (t-1)*wⱼᶜ+Lⱼ ≤ t)
        end

        # Sᵁ_dash = Σ(μₖ*wᵢᶜ+lᵢ)
        # Sᴸ_dash = Σ(μₖ*wᵢᶜ-lᵢ)
        Sᵁ_dash = sum(map(i -> (t-1)*Wᶜ[i,k]+L[i], filter(i -> i != k, 1:n)))
        Sᴸ_dash = sum(map(i -> (t-1)*Wᶜ[i,k]-L[i], filter(i -> i != k, 1:n)))
        @constraint(model, Sᵁ_dash + 1-Lₖ ≥ t)
        @constraint(model, Sᴸ_dash + 1+Lₖ ≤ t)
        @constraint(model, 1-Lₖ ≥ t*ε)
        Σl = sum(map(j-> L[j], filter(j -> j != k, 1:n)))
        @constraint(model, Σl ≤ (t-1)*d⃰+ε)

        # dₖ  = sum(map(j -> l[j], filter(j -> j != k, 1:n)))
        @objective(model, Min, L[k])

        optimize!(model)
        # println("k=$k, t = ", value.(t))
        tₖ⃰ = value.(t)
        L⃰ = value.(L) 
        return (
            tₖ⃰ = tₖ⃰ ,
            L⃰ = L⃰
        )

    finally
        empty!(model)
    end
end

# 提案手法 AMR-E, AMR-G, AMR-A
@inline function AMR_kai(A::Matrix{T}, method::Function)::AMRkai_Individual{T} where {T <: Real}

    if !isCrispPCM(A)
        throw(ArgumentError("A is not a crisp PCM"))
    end

    m, n = size(A)

    # Phase 1
    Wᶜ = AMRkai_phase1(A, method)

    # Phase 2
    d⃰ = Vector{T}(undef, n) 
    for k = 1:n
        d⃰[k] = AMRkai_phase2_jump(A, Wᶜ, k, n)
    end

    # Phase 3
    wᴸ = Matrix{T}(undef, m, n)
    wᵁ = Matrix{T}(undef, m, n)
    centers = Matrix{T}(undef, m, n)
    l = Matrix{T}(undef, m, n)

    for k = 1:n
        (tₖ⃰, L⃰) = AMRkai_phase3_jump(A, Wᶜ, d⃰[k], k, n)

        for i = 1:n
            lᵢ⃰ = L⃰[i]/tₖ⃰
            wᵢᶜ = Wᶜ[i,k]
            if i==k
                wᴸᵢ = 1/tₖ⃰ - lᵢ⃰
                wᵁᵢ = 1/tₖ⃰ + lᵢ⃰
                centers[i, k] = 1/tₖ⃰
                l[i, k] = lᵢ⃰

            else
                wᴸᵢ = (1-1/tₖ⃰)*wᵢᶜ - lᵢ⃰
                wᵁᵢ = (1-1/tₖ⃰)*wᵢᶜ + lᵢ⃰
                # println(tₖ⃰)
                centers[i, k] = (1-1/tₖ⃰)*wᵢᶜ
                # wᴸᵢ = wᵢᶜ - lᵢ⃰
                # wᵁᵢ = wᵢᶜ + lᵢ⃰
                # centers[i, k] = wᵢᶜ
                l[i, k] = lᵢ⃰
            end
            # 各kの時の推定値を縦ベクトルとして格納
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

    w_c = sum(w̅̅ᴸ.+ w̅̅ᵁ)/2


    for i = 1:n
        W̅̅[i] = (w̅̅ᴸ[i])..(w̅̅ᵁ[i])
    end


    return (
        s=w_c,
        centers=centers, l=l,
        wᴸ=w̅̅ᴸ, wᵁ=w̅̅ᵁ,
        W=W̅̅
    )

end