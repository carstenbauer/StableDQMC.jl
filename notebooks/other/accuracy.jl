cd(@__DIR__); using Pkg; Pkg.activate(".");

using Revise, StableDQMC, GenericSVD, JacobiSVD, LinearAlgebra, Plots, Statistics, LaTeXStrings, DataStructures
using BenchmarkTools, StatsPlots, DataFrames, Infiltrator

pyplot();
# PyPlot.plt.style.use("default")
PyPlot.plt.style.use("publication_tex")
default(
    framestyle = :grid,
    grid = true,
    guidefont=font(15),
#     framestyle = :box,
#     grid = false,
#     guidefont=font(18),
    legend = true,
    size = (1.6*360, 360),
    xtickfont=font(13),
    ytickfont=font(13),
    legendfont=font(12),
    linewidth = 1.5
);

## Model

# Hubbard
L = 2
nsites = L^2
n = nsites * 2 # spin
t = 1
U = 6
μ = 0
Δτ = 0.1

# Kinetic part (PBC)
T = diagm(-1 => fill(-t, n-1), 1 => fill(-t, n-1), 0 => fill(μ, n))
T[1,end] = -t; T[end,1] = -t;

# on-site interaction
λ = acosh(exp(0.5*U*Δτ))
eV = Diagonal(vcat(fill(exp(λ), nsites), fill(exp(-λ), nsites)))

B = exp(-Δτ .* T) * eV

println()
@info "Hubbard model: U = $U, L = $L"
@info "Condition number of B: $(cond(B))"

## Bchain product

Bbig = BigFloat.(B)
N = 500

setprecision(100) do
    global svs_genericsvd = calc_Bchain_svd(Bbig, N; svdalg = genericsvd)[2];
end
svs_qr_udt = calc_Bchain_qr(B, N)[2];
svs_gesvd  = calc_Bchain_svd(B, N; svdalg = gesvd)[2];
svs_gesdd  = calc_Bchain_svd(B, N; svdalg = gesdd)[2];
svs_gesvj  = calc_Bchain_svd(B, N; svdalg = gesvj)[2];

betas = range(1, N * 0.1, length=N)
plot(betas, svs_genericsvd, color="#CA0020", labels = reshape(["exact", "", "","", "","", "", ""], (1,8)))
plot!(betas, svs_qr_udt, color="#E66101", labels = reshape(["QR", "", "","", "","", "", ""], (1,8)))
plot!(betas, svs_gesvj, color="#FDB863", labels = reshape(["SVD (Jacobi)", "", "","", "","", "", ""], (1,8)))
plot!(betas, svs_gesvd, color="#008837", labels = reshape(["SVD", "", "","", "","", "", ""], (1,8)))
plot!(betas, svs_gesdd, color="#5E3C99", labels = reshape(["SVD (D\\&C)", "", "","", "","", "", ""], (1,8)))

ylabel!("log singular values")
xlabel!(L"inverse temperature $\beta$")
# ylims!(-160, 100)
xlims!(0., 50.)
# savefig("../paper/figures/decomp_comparison_simple.pdf")

## Equal-time Green's function

function accuracy_greens(B; svdinversion = inv_one_plus)
    Bbig = BigFloat.(B)

    Δs_qr    = SortedDict{Int, Float64}()
    Δs_gesvd = SortedDict{Int, Float64}()
    Δs_gesdd = SortedDict{Int, Float64}()
    Δs_gesvj = SortedDict{Int, Float64}()

    for beta in range(5, 40, step=5)
        N = Int(beta * 10)

        local G_genericsvd
        setprecision(1000) do
            F_genericsvd = calc_Bchain_svd(Bbig, N; svdalg = genericsvd)[1];
            G_genericsvd = inv(I + Matrix(F_genericsvd))
        end

        # Bchain
        F_qr_udt = calc_Bchain_qr(B, N)[1];
        F_gesvd  = calc_Bchain_svd(B, N; svdalg = gesvd)[1];
        F_gesdd  = calc_Bchain_svd(B, N; svdalg = gesdd)[1];
        F_gesvj  = calc_Bchain_svd(B, N; svdalg = StableDQMC.gesvj)[1];

        # inversion
        G_qr_udt = inv_one_plus(F_qr_udt)
        G_gesvd  = svdinversion(F_gesvd; svdalg = gesvd!)
        G_gesdd  = svdinversion(F_gesdd; svdalg = gesdd!)
        G_gesvj  = svdinversion(F_gesvj; svdalg = StableDQMC.gesvj!)

        Δs_qr[beta]    = maximum(abs.(G_qr_udt - G_genericsvd))
        Δs_gesvd[beta] = maximum(abs.(G_gesvd - G_genericsvd))
        Δs_gesdd[beta] = maximum(abs.(G_gesdd - G_genericsvd))
        Δs_gesvj[beta] = maximum(abs.(G_gesvj - G_genericsvd))
    end


    return Δs_qr, Δs_gesvd, Δs_gesdd, Δs_gesvj
end

# Δs_qr, Δs_gesvd, Δs_gesdd, Δs_gesvj = Juno.@enter accuracy_greens(B);
Δs_qr, Δs_gesvd, Δs_gesdd, Δs_gesvj = accuracy_greens(B);

##

betas = collect(keys(Δs_qr))
p = plot(betas, log.(10, values(Δs_qr)), color="#E66101", marker=true, label="QR")
plot!(betas, log.(10, values(Δs_gesvd)), color="#008837", marker=true, label="SVD")
plot!(betas, log.(10, values(Δs_gesdd)), color="#5E3C99", marker=true, label="SVD (D\\&C)")
plot!(betas, log.(10, values(Δs_gesvj)), color="#FDB863", marker=true, label="SVD (Jacobi)")

xlabel!(L"inverse temperature $\beta$")
ylabel!(L"\log(\textrm{max}(\textrm{abs}(G - G_{\textrm{exact}})))")
# ylims!(0., 1e-3)
savefig("../paper/figures/accuracy_greens_svd_regularinv_U_$(U)_L_$(L)_t_$(t).pdf")

## Time-displaced Green's function

# G(slice, 1) = [B(slice, 1)^-1 + B(beta, slice)]^-1
# G(slice, 1) = [B(slice, 1)^-1 + B(beta, slice)]^-1
function accuracy_tdgf(B; loh = false)
    Bbig = BigFloat.(B)

#     betas = range(5, 40, step=5)
    beta = 40
    M = Int(beta * 10)
    slices = 1:10:M

    df = DataFrame(slice=Int[], qr=Float64[], gesdd=Float64[], gesvd=Float64[], gesvj=Float64[])

    inv_sum_method = loh ? inv_sum_loh : inv_sum

    for s in slices
        local TDGF_genericsvd
        setprecision(1000) do
            TDGF_genericsvd = calc_tdgf_svd(Bbig, s, M-s, inv_sum_method = inv_sum_method, svdalg_chain = genericsvd, svdalg_inv = genericsvd!)
        end

        TDGF_qr = calc_tdgf_qr(B, s, M-s, inv_sum_method = inv_sum_method)
        TDGF_gesdd = calc_tdgf_svd(B, s, M-s, inv_sum_method = inv_sum_method, svdalg_chain = gesdd, svdalg_inv = gesdd!)
        TDGF_gesvd = calc_tdgf_svd(B, s, M-s, inv_sum_method = inv_sum_method, svdalg_chain = gesvd, svdalg_inv = gesvd!)
        TDGF_gesvj = calc_tdgf_svd(B, s, M-s, inv_sum_method = inv_sum_method, svdalg_chain = StableDQMC.gesvj, svdalg_inv = StableDQMC.gesvj!)

        Δ_qr    = maximum(abs, TDGF_qr - TDGF_genericsvd)
        Δ_gesdd = maximum(abs, TDGF_gesdd - TDGF_genericsvd)
        Δ_gesvd = maximum(abs, TDGF_gesvd - TDGF_genericsvd)
        Δ_gesvj = maximum(abs, TDGF_gesvj - TDGF_genericsvd)

        push!(df, [s, Δ_qr, Δ_gesdd, Δ_gesvd, Δ_gesvj])
    end

    return df
end

df_tdgf = accuracy_tdgf(B; loh = true)

##

@df df_tdgf plot(:slice .* 0.1, (x -> log.(10, x)).([:qr, :gesvd, :gesdd, :gesvj]),
            marker=true,
            color=permutedims(["#E66101", "#008837", "#5E3C99", "#FDB863"]), # #E66101 for QR
            label=permutedims(["QR", "SVD", "SVD (D\\&C)", "SVD (Jacobi)"]),
#             xlims=(0,10)
           )
xlabel!(L"imaginary time $\tau$")
ylabel!(L"\log(\textrm{max}(\textrm{abs}(G(\tau, 0) - G_{\textrm{exact}}(\tau, 0))))")
savefig("../paper/figures/accuracy_tdgf_loh_U_$(U)_L_$(L)_t_$(t).pdf")
