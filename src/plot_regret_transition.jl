using Plots
using LinearAlgebra

function plot_regret_transition(utility::Matrix{Float64}, methodW::Vector{NamedTuple{(:L, :R, :adjacent), Tuple{Vector{Float64}, Vector{Float64}, Float64}}}, t_range::Tuple{Float64, Float64})
    # リグレット行列の初期化
    matrix = create_minimax_R_Matrix(utility)
    
    # t値の範囲を生成（より細かい分割）
    t_values = range(t_range[1], t_range[2], length=200)
    
    # 各代替案のリグレットを格納する配列
    regrets = zeros(5, length(t_values))
    
    # 各t値でリグレットを計算
    for (i, t) in enumerate(t_values)
        Y_L = methodW[1].L .* t
        Y_R = methodW[1].R .* t
        
        regret_max, _ = calc_regret(matrix, Y_L, Y_R)
        regrets[:, i] = [r[1] for r in regret_max]
    end
    
    # プロット作成
    p = plot(
        t_values,
        regrets',
        xlabel="t",
        ylabel="Regret",
        title="Regret Transition",
        lw=2,
        legend=:outerright,  # 凡例を右側の外に配置
        size=(800,400),      # プロットのサイズを調整して凡例用のスペースを確保
        grid=true
    )
    
    # t_rangeの範囲を点線で表示
    vline!([t_range[1], t_range[2]], 
           label="t range", 
           linestyle=:dash, 
           color=:black, 
           alpha=0.5)
    
    return p
end
