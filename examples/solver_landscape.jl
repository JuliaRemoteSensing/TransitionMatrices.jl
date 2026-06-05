### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 246053ac-6082-11f1-bed7-11580218e966
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(; path = dirname(@__DIR__))
    Pkg.instantiate()
    using TransitionMatrices, Plots
end

# ╔═╡ 2460532a-6082-11f1-83c5-3796c2e835e2
md"""
# Solver landscape & convergence

The same T-matrix can be built several ways. The two-layer API selects them with
a solver object passed to `transition_matrix(s, λ, solver)`:

- a **fixed** solver — `EBCM(nₘₐₓ, Ng)`, `IITM(nₘₐₓ, Nr, Nϑ)` — one build;
- an **iterative** solver — `Iterative(EBCM; …)` — sweeps to convergence;

plus the `stable = true` option that removes the high-aspect EBCM cancellation.
This notebook shows *why* it matters: on a high-aspect spheroid, classic EBCM
converges briefly and then **diverges** as `nₘₐₓ` grows, while the stabilized
path stays converged.
"""

# ╔═╡ 246053b6-6082-11f1-a95d-e7ab43aaff16
md"""
## A high-aspect prolate spheroid (aspect ≈ 4)

Three routes to the same cross section:
"""

# ╔═╡ 246053c2-6082-11f1-9425-cd50e8e63bdd
begin
    λ = 2π
    prolate = TransitionMatrices.Spheroid{Float64, ComplexF64}(2.5198421, 10.079368, 1.55 + 0.01im)

    routes = [
        (; route = "EBCM fixed (n=24)", Csca = calc_Csca(transition_matrix(prolate, λ, EBCM(24, 384)), λ)),
        (; route = "Iterative(EBCM; stable)", Csca = calc_Csca(transition_matrix(prolate, λ, Iterative(EBCM; stable = true)), λ)),
        (; route = "IITM", Csca = calc_Csca(transition_matrix(prolate, λ, IITM(24, 30, 60)), λ)),
    ]
end

# ╔═╡ 246053ca-6082-11f1-ba56-6f4f921058ef
md"""
## Convergence vs truncation order

Relative error of ``Q_\text{sca}`` against a high-order stabilized reference, for
classic EBCM vs the `stable = true` EBCM. Watch the classic curve blow up.
"""

# ╔═╡ 246053d4-6082-11f1-b47b-0143b9abc039
begin
    ref = calc_Csca(transition_matrix(prolate, λ, EBCM(40, 600; stable = true)), λ)
    ns = 16:2:40
    classic_err = [abs(calc_Csca(transition_matrix(prolate, λ, EBCM(n, 16n)), λ) - ref) / ref for n in ns]
    stable_err = [abs(calc_Csca(transition_matrix(prolate, λ, EBCM(n, 16n; stable = true)), λ) - ref) / ref for n in ns]
    nothing
end

# ╔═╡ 246053de-6082-11f1-a6ce-e998a04e4ad7
let
    plot(collect(ns), classic_err; yscale = :log10, lw = 2, marker = :circle,
        label = "classic EBCM", legend = :left)
    plot!(collect(ns), stable_err; lw = 2, marker = :square, label = "stable EBCM")
    plot!(xlabel = "truncation order nₘₐₓ", ylabel = "relative error in Qsca",
        title = "high-aspect spheroid: classic EBCM diverges, stable holds",
        size = (720, 420))
end

# ╔═╡ 246053fc-6082-11f1-9730-eb288864dde4
md"""
## Choosing a solver

| situation | solver |
|---|---|
| axisymmetric, moderate aspect | `calc_T(s, λ)` (auto-converged classic EBCM) |
| high-aspect **spheroid** | `Iterative(EBCM; stable = true)` |
| cylinder / Chebyshev / 3-D, hard cases | `IITM(nₘₐₓ, Nr, Nϑ[, Nφ])` |
| wavelength / index **sweep** | `ShMatrix` / `prepare_sh` |
| explicit discretization (benchmarking) | `EBCM(nₘₐₓ, Ng)`, `IITM(…)` |

The iterative layer (`Iterative`) wraps any of these with automatic convergence.
"""

# ╔═╡ Cell order:
# ╟─2460532a-6082-11f1-83c5-3796c2e835e2
# ╠═246053ac-6082-11f1-bed7-11580218e966
# ╟─246053b6-6082-11f1-a95d-e7ab43aaff16
# ╠═246053c2-6082-11f1-9425-cd50e8e63bdd
# ╟─246053ca-6082-11f1-ba56-6f4f921058ef
# ╠═246053d4-6082-11f1-b47b-0143b9abc039
# ╠═246053de-6082-11f1-a6ce-e998a04e4ad7
# ╟─246053fc-6082-11f1-9730-eb288864dde4
