using IntervalArithmetic
using JuMP
import HiGHS

include("./crisp-pcm.jl")
include("./nearly-equal.jl")
include("./solve-deterministic-ahp.jl")

"""
従来法のMMR-Wをjuliaにて実装
提案手法の実装が正しいことを確認するために、先行研究の結果と比較する
"""


LPResult_Individual = @NamedTuple{
    # 区間重みベクトル
    wᴸ::Vector{T}, wᵁ::Vector{T},
    W::Vector{Interval{T}}, # ([wᵢᴸ, wᵢᵁ])
} where {T<:Real}

# Phase2のループの中の部分
@inline function phase2_jump(A::Matrix{T}, Wᶜ::Vector{T}, k::Int, n::Int)::T where {T<:Real}
    ε = 1e-6 # << 1

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, l[i=1:n] ≥ ε)

        for j = 1:n
            for i = filter(i -> i != j, 1:n)
                aᵢⱼ = A[i, j]
                wᵢᶜ = Wᶜ[i]
                wⱼᶜ = Wᶜ[j]
                lᵢ = l[i]
                lⱼ = l[j]
                @constraint(model, aᵢⱼ * (wⱼᶜ - lⱼ) ≤ wᵢᶜ + lᵢ)
            end
        end

        for j = 1:n
            lⱼ = l[j]
            wⱼᶜ = Wᶜ[j]
            # Sᵁ = Σ(μₖ*wᵢᶜ+lᵢ)
            # Sᴸ = Σ(μₖ*wᵢᶜ-lᵢ)
            Sᵁ = sum(map(i -> Wᶜ[i] + l[i], filter(i -> i != j, 1:n)))
            Sᴸ = sum(map(i -> Wᶜ[i] - l[i], filter(i -> i != j, 1:n)))
            @constraint(model, Sᵁ + wⱼᶜ - lⱼ ≥ 1)
            @constraint(model, Sᴸ + wⱼᶜ + lⱼ ≤ 1)
        end

        for i = 1:n
            wᵢᶜ = Wᶜ[i]
            lᵢ = l[i]
            @constraint(model, wᵢᶜ - lᵢ ≥ ε)
        end

        dₖ = sum(map(j -> l[j], filter(j -> j != k, 1:n)))
        @objective(model, Min, dₖ)

        optimize!(model)
        d̂ₖ = sum(map(j -> value.(l[j]), filter(j -> j != k, 1:n)))

        return d̂ₖ

    finally
        empty!(model)
    end
end

# Phase3のループの中の部分
@inline function phase3_jump(A::Matrix{T}, Wᶜ::Vector{T}, d̂::T, k::Int, n::Int)::Vector{T} where {T<:Real}
    ε = 1e-6 # << 1

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, l[i=1:n] ≥ ε)
        lₖ = l[k]

        for j = 1:n
            for i = filter(i -> i != j, 1:n)
                aᵢⱼ = A[i, j]
                wᵢᶜ = Wᶜ[i]
                wⱼᶜ = Wᶜ[j]
                lᵢ = l[i]
                lⱼ = l[j]
                @constraint(model, aᵢⱼ * (wⱼᶜ - lⱼ) ≤ wᵢᶜ + lᵢ)
            end
        end

        for j = 1:n
            lⱼ = l[j]
            wⱼᶜ = Wᶜ[j]
            # Sᵁ = Σ(μₖ*wᵢᶜ+lᵢ)
            # Sᴸ = Σ(μₖ*wᵢᶜ-lᵢ)
            Sᵁ = sum(map(i -> Wᶜ[i] + l[i], filter(i -> i != j, 1:n)))
            Sᴸ = sum(map(i -> Wᶜ[i] - l[i], filter(i -> i != j, 1:n)))
            @constraint(model, Sᵁ + wⱼᶜ - lⱼ ≥ 1)
            @constraint(model, Sᴸ + wⱼᶜ + lⱼ ≤ 1)
        end

        for i = 1:n
            wᵢᶜ = Wᶜ[i]
            lᵢ = l[i]
            @constraint(model, wᵢᶜ - lᵢ ≥ ε)
        end

        Σl = sum(map(j -> l[j], filter(j -> j != k, 1:n)))
        @constraint(model, Σl ≤ d̂ + ε)

        @objective(model, Min, lₖ)

        optimize!(model)
        l̂ = value.(l)
        return l̂

    finally
        empty!(model)
    end
end

# 提案手法
@inline function MMR_W(A::Matrix{T}, method::Function)::LPResult_Individual{T} where {T<:Real}

    if !isCrispPCM(A)
        throw(ArgumentError("A is not a crisp PCM"))
    end

    m, n = size(A)

    # Phase 1
    # EV, GMで区間重要度の中心を求める
    Wᶜ = method(A)

    # Phase 2
    d̂ = Vector{T}(undef, n)
    for k = 1:n
        d̂[k] = phase2_jump(A, Wᶜ, k, n)
    end

    # Phase 3
    wᴸ = Matrix{T}(undef, m, n)
    wᵁ = Matrix{T}(undef, m, n)
    l̂ = Matrix{T}(undef, m, n)
    for k = 1:n
        l = phase3_jump(A, Wᶜ, d̂[k], k, n)
        for j = 1:n
            l̂[j, k] = l[j]
        end
    end

    # Phase 4
    w̅̅ᴸ = Vector{T}(undef, n)
    w̅̅ᵁ = Vector{T}(undef, n)
    W̅̅ = Vector{Interval{T}}(undef, n)

    for i = 1:n
        wᵢᶜ = Wᶜ[i]
        w̅̅ᴸ[i] = wᵢᶜ - maximum(l̂[i, :])
        w̅̅ᵁ[i] = wᵢᶜ + maximum(l̂[i, :])

        # precision error 対応
        if w̅̅ᴸ[i] > w̅̅ᵁ[i]
            w̅̅ᴸ[i] = w̅̅ᵁ[i]
        end

        W̅̅[i] = interval(w̅̅ᴸ[i], w̅̅ᵁ[i])
    end

    return (
        wᴸ=w̅̅ᴸ, wᵁ=w̅̅ᵁ,
        W=W̅̅
    )

end

#-----------------------------------------------------------------------
# ここからW/C法の実装（追記部分）
#-----------------------------------------------------------------------

"""
Phase2 (W/C法バージョン)
目的関数を Σ(lᵢ/wᵢᶜ) に変更
"""
@inline function phase2_jump_wc(A::Matrix{T}, Wᶜ::Vector{T}, k::Int, n::Int)::T where {T<:Real}
    ε = 1e-6 # << 1

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, l[i=1:n] ≥ ε)

        # 制約条件はMMR_Wと同じ
        for j = 1:n
            for i = filter(i -> i != j, 1:n)
                aᵢⱼ = A[i, j]
                wᵢᶜ = Wᶜ[i]
                wⱼᶜ = Wᶜ[j]
                lᵢ = l[i]
                lⱼ = l[j]
                @constraint(model, aᵢⱼ * (wⱼᶜ - lⱼ) ≤ wᵢᶜ + lᵢ)
            end
        end

        for j = 1:n
            lⱼ = l[j]
            wⱼᶜ = Wᶜ[j]
            Sᵁ = sum(map(i -> Wᶜ[i] + l[i], filter(i -> i != j, 1:n)))
            Sᴸ = sum(map(i -> Wᶜ[i] - l[i], filter(i -> i != j, 1:n)))
            @constraint(model, Sᵁ + wⱼᶜ - lⱼ ≥ 1)
            @constraint(model, Sᴸ + wⱼᶜ + lⱼ ≤ 1)
        end

        for i = 1:n
            wᵢᶜ = Wᶜ[i]
            lᵢ = l[i]
            @constraint(model, wᵢᶜ - lᵢ ≥ ε)
        end

        # ★★★ 目的関数を Σ(l/wᶜ) に変更 ★★★
        dₖ_wc = sum(map(j -> l[j] / Wᶜ[j], filter(j -> j != k, 1:n)))
        @objective(model, Min, dₖ_wc)

        optimize!(model)
        # 最適値も目的関数に合わせて計算
        d̂ₖ_wc = sum(map(j -> value.(l[j]) / Wᶜ[j], filter(j -> j != k, 1:n)))

        return d̂ₖ_wc

    finally
        empty!(model)
    end
end

"""
Phase3 (W/C法バージョン)
目的関数を lₖ/wᶜₖ に変更
"""
@inline function phase3_jump_wc(A::Matrix{T}, Wᶜ::Vector{T}, d̂_wc::T, k::Int, n::Int)::Vector{T} where {T<:Real}
    ε = 1e-6 # << 1

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    try
        @variable(model, l[i=1:n] ≥ ε)
        lₖ = l[k]

        # 制約条件はMMR_Wと同じ
        for j = 1:n
            for i = filter(i -> i != j, 1:n)
                aᵢⱼ = A[i, j]
                wᵢᶜ = Wᶜ[i]
                wⱼᶜ = Wᶜ[j]
                lᵢ = l[i]
                lⱼ = l[j]
                @constraint(model, aᵢⱼ * (wⱼᶜ - lⱼ) ≤ wᵢᶜ + lᵢ)
            end
        end

        for j = 1:n
            lⱼ = l[j]
            wⱼᶜ = Wᶜ[j]
            Sᵁ = sum(map(i -> Wᶜ[i] + l[i], filter(i -> i != j, 1:n)))
            Sᴸ = sum(map(i -> Wᶜ[i] - l[i], filter(i -> i != j, 1:n)))
            @constraint(model, Sᵁ + wⱼᶜ - lⱼ ≥ 1)
            @constraint(model, Sᴸ + wⱼᶜ + lⱼ ≤ 1)
        end

        for i = 1:n
            wᵢᶜ = Wᶜ[i]
            lᵢ = l[i]
            @constraint(model, wᵢᶜ - lᵢ ≥ ε)
        end

        # ★★★ Phase2の最適値に関する制約も Σ(l/wᶜ) で評価 ★★★
        Σl_wc = sum(map(j -> l[j] / Wᶜ[j], filter(j -> j != k, 1:n)))
        @constraint(model, Σl_wc ≤ d̂_wc + ε)

        # ★★★ 目的関数を lₖ / wᶜₖ に変更 ★★★
        @objective(model, Min, lₖ / Wᶜ[k])

        optimize!(model)
        l̂ = value.(l)
        return l̂

    finally
        empty!(model)
    end
end


"""
MMR-W/C法
"""
@inline function MMR_WC(A::Matrix{T}, method::Function)::LPResult_Individual{T} where {T<:Real}

    if !isCrispPCM(A)
        throw(ArgumentError("A is not a crisp PCM"))
    end

    m, n = size(A)

    # Phase 1: 区間中心の計算（変更なし）
    Wᶜ = method(A)

    # Phase 2: d̂の計算（w/c版の関数を呼ぶ）
    d̂_wc = Vector{T}(undef, n)
    for k = 1:n
        d̂_wc[k] = phase2_jump_wc(A, Wᶜ, k, n)
    end

    # Phase 3: 最適なlをn通り計算（w/c版の関数を呼ぶ）
    l̂ = Matrix{T}(undef, m, n)
    for k = 1:n
        l = phase3_jump_wc(A, Wᶜ, d̂_wc[k], k, n)
        for j = 1:n
            l̂[j, k] = l[j]
        end
    end

    # Phase 4: MMRによる統合（変更なし）
    w̅̅ᴸ = Vector{T}(undef, n)
    w̅̅ᵁ = Vector{T}(undef, n)
    W̅̅ = Vector{Interval{T}}(undef, n)

    for i = 1:n
        wᵢᶜ = Wᶜ[i]
        w̅̅ᴸ[i] = wᵢᶜ - maximum(l̂[i, :])
        w̅̅ᵁ[i] = wᵢᶜ + maximum(l̂[i, :])

        if w̅̅ᴸ[i] > w̅̅ᵁ[i]
            w̅̅ᴸ[i] = w̅̅ᵁ[i]
        end

        W̅̅[i] = interval(w̅̅ᴸ[i], w̅̅ᵁ[i])
    end

    return (wᴸ=w̅̅ᴸ, wᵁ=w̅̅ᵁ, W=W̅̅)
end

"""
AMR-W/C法
"""
@inline function AMR_WC(A::Matrix{T}, method::Function)::LPResult_Individual{T} where {T<:Real}

    if !isCrispPCM(A)
        throw(ArgumentError("A is not a crisp PCM"))
    end

    m, n = size(A)

    # Phase 1: 区間中心の計算（変更なし）
    Wᶜ = method(A)

    # Phase 2: d̂の計算（w/c版の関数を呼ぶ）
    d̂_wc = Vector{T}(undef, n)
    for k = 1:n
        d̂_wc[k] = phase2_jump_wc(A, Wᶜ, k, n)
    end

    # Phase 3: 最適なlをn通り計算（w/c版の関数を呼ぶ）
    l̂ = Matrix{T}(undef, m, n)
    for k = 1:n
        l = phase3_jump_wc(A, Wᶜ, d̂_wc[k], k, n)
        for j = 1:n
            l̂[j, k] = l[j]
        end
    end

    # Phase 4: ★★★ AMRによる統合（meanを使用） ★★★
    w̅̅ᴸ = Vector{T}(undef, n)
    w̅̅ᵁ = Vector{T}(undef, n)
    W̅̅ = Vector{Interval{T}}(undef, n)

    for i = 1:n
        wᵢᶜ = Wᶜ[i]
        # maximum の代わりに sum/n (mean) を使う
        avg_l = sum(l̂[i, :]) / n
        w̅̅ᴸ[i] = wᵢᶜ - avg_l
        w̅̅ᵁ[i] = wᵢᶜ + avg_l

        if w̅̅ᴸ[i] > w̅̅ᵁ[i]
            w̅̅ᴸ[i] = w̅̅ᵁ[i]
        end

        W̅̅[i] = interval(w̅̅ᴸ[i], w̅̅ᵁ[i])
    end

    return (wᴸ=w̅̅ᴸ, wᵁ=w̅̅ᵁ, W=W̅̅)
end