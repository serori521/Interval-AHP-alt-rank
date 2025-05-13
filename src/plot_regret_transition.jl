function plot_regret_transition(
    utility::Matrix{Float64},
    methodW_vector::Vector{<:NamedTuple},
    t_range::Tuple{Float64,Float64}
)
    # リグレット行列の初期化
    matrix = create_minimax_R_Matrix(utility)

    # t値の範囲を生成
    t_values = range(t_range[1], t_range[2], length=200)

    n_alternatives = size(utility, 1)

    # 各代替案の最大リグレットを格納する配列
    regrets_over_t = zeros(n_alternatives, length(t_values))

    if isempty(methodW_vector)
        println("警告: plot_regret_transition - methodW_vector が空です。")
        return plot(title="各代替案の最大リグレット推移 (データ不足)", xlabel="t", ylabel="最大リグレット値", size=(800, 500))
    end
    current_methodW = methodW_vector[1] # Vectorの最初の要素を使用

    # 各t値でリグレットを計算
    for (col_idx, t) in enumerate(t_values)
        Y_L = current_methodW.L .* t
        Y_R = current_methodW.R .* t

        temp_matrix_for_calc = deepcopy(matrix)
        regret_max_output, _ = calc_regret(temp_matrix_for_calc, Y_L, Y_R)

        if length(regret_max_output) == n_alternatives
            try
                # 元のコードの [r[1] for r in regret_max] のパターンを想定
                regrets_over_t[:, col_idx] = [r[1] for r in regret_max_output]
            catch e
                if eltype(regret_max_output) <: Real # 数値のベクトルの場合
                    regrets_over_t[:, col_idx] = regret_max_output
                else # それ以外はエラーとしてNaN
                    if col_idx == 1
                        println("警告(plot_regret_transition): t=$t での calc_regret の戻り値形式が不正。詳細: $e")
                    end
                    regrets_over_t[:, col_idx] .= NaN
                end
            end
        else
            if col_idx == 1
                println("警告(plot_regret_transition): t=$t での calc_regret の戻り値の要素数が不一致。")
            end
            regrets_over_t[:, col_idx] .= NaN
        end
    end

    # プロット作成
    p = plot(
        xlabel="t",
        ylabel="最大リグレット値",
        title="各代替案の最大リグレット推移",
        legend=:outerright,
        size=(800, 500),
        grid=true
    )

    colors = [:red, :blue, :green, :purple, :orange, :brown, :pink, :gray, :cyan, :magenta]

    for alt_idx in 1:n_alternatives
        plot!(p, t_values, regrets_over_t[alt_idx, :],
            label="代替案 $alt_idx",
            color=colors[alt_idx > length(colors) ? mod1(alt_idx, length(colors)) : alt_idx],
            linewidth=2
        )
    end

    vline!(p, [t_range[1], t_range[2]],
        label="t の分析範囲",
        linestyle=:dash,
        color=:black,
        alpha=0.5
    )

    return p
end