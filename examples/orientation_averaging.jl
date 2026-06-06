### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 2366dfd4-6082-11f1-bc99-e3b4423633f0
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(; path = dirname(@__DIR__))
    Pkg.instantiate()
    using TransitionMatrices, Plots
end

# ╔═╡ 2366df48-6082-11f1-880f-99c1a49d5d94
md"""
# Orientation averaging

For randomly oriented particles you need the orientation-averaged scattering
quantities. `TransitionMatrices.jl` offers two routes:

- **`RandomOrientationTransitionMatrix(T)`** — Mishchenko's *analytic* closed
  form for the orientation average (fast, exact);
- **`orientation_average(T, p; Nα, Nβ, Nγ)`** — *numerical* average over the
  Euler angles against an orientation distribution `p`.

This notebook checks they agree and shows the numerical average converging to the
analytic value as the angular grid is refined.
"""

# ╔═╡ 2366dfe8-6082-11f1-b905-b7bae29b9c75
md"""
## A non-spherical particle
"""

# ╔═╡ 2366dff2-6082-11f1-bcad-09c1d7464bc5
begin
    λ = 2π
    spheroid = TransitionMatrices.Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)
    T = calc_T(spheroid, λ)
    uniform = (α, β, γ) -> 1 / (8π^2)   # isotropic orientation distribution
end

# ╔═╡ 2366dff2-6082-11f1-b19f-9f43b5fe61a4
md"""
## Analytic vs numerical

Both should give the same orientation-averaged efficiencies.
"""

# ╔═╡ 2366dffc-6082-11f1-8c98-1fba91e58dc3
begin
    Tanalytic = RandomOrientationTransitionMatrix(T)
    Tnumeric = orientation_average(T, uniform; Nα = 20, Nβ = 20, Nγ = 20)

    compare = [
        (; method = "analytic (Mishchenko)", Qsca = calc_Csca(Tanalytic, λ),
            Qext = calc_Cext(Tanalytic, λ), g = asymmetry_parameter(Tanalytic, λ)),
        (; method = "numerical (20³ grid)", Qsca = calc_Csca(Tnumeric, λ),
            Qext = calc_Cext(Tnumeric, λ), g = asymmetry_parameter(Tnumeric, λ))
    ]
end

# ╔═╡ 2366e008-6082-11f1-b755-1527ca352dc3
md"""
## Convergence of the numerical average

The numerical average approaches the analytic value as the `β`-grid is refined;
the analytic closed form gives it directly, at a fraction of the cost.
"""

# ╔═╡ 2366e024-6082-11f1-bc7a-5d2a795a9fec
begin
    Qref = calc_Csca(Tanalytic, λ)
    Ns = 4:2:18
    errs = [abs(calc_Csca(orientation_average(T, uniform; Nα = 2, Nβ = N, Nγ = 2), λ) -
                Qref) / Qref
            for N in Ns]
    nothing
end

# ╔═╡ 2366e02e-6082-11f1-9678-c1ad121a4222
let
    plot(collect(Ns), errs; yscale = :log10, lw = 2, marker = :circle,
        xlabel = "β-grid points Nβ", ylabel = "|Qsca_numeric − Qsca_analytic| / Qsca",
        label = "numerical average", legend = :topright,
        title = "convergence to the analytic orientation average", size = (680, 400))
end

# ╔═╡ 2366e03a-6082-11f1-958b-71200c5d7a31
md"""
## Takeaway

The analytic `RandomOrientationTransitionMatrix` is both exact and much cheaper
than refining a 3-D Euler-angle quadrature — prefer it for random-orientation
cross sections. The numerical `orientation_average` remains useful for
*non-uniform* orientation distributions (just pass a different `p`).
"""

# ╔═╡ Cell order:
# ╟─2366df48-6082-11f1-880f-99c1a49d5d94
# ╠═2366dfd4-6082-11f1-bc99-e3b4423633f0
# ╟─2366dfe8-6082-11f1-b905-b7bae29b9c75
# ╠═2366dff2-6082-11f1-bcad-09c1d7464bc5
# ╟─2366dff2-6082-11f1-b19f-9f43b5fe61a4
# ╠═2366dffc-6082-11f1-8c98-1fba91e58dc3
# ╟─2366e008-6082-11f1-b755-1527ca352dc3
# ╠═2366e024-6082-11f1-bc7a-5d2a795a9fec
# ╠═2366e02e-6082-11f1-9678-c1ad121a4222
# ╟─2366e03a-6082-11f1-958b-71200c5d7a31
