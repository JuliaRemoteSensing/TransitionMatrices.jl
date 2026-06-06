### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 162c8ff8-6080-11f1-af6b-f5508490c007
begin
    import Pkg
    Pkg.activate(@__DIR__)                 # the examples/ environment
    Pkg.develop(; path = dirname(@__DIR__)) # point TransitionMatrices at this repo
    Pkg.instantiate()
    using TransitionMatrices, ForwardDiff, Plots, BenchmarkTools
end

# ╔═╡ 162c8f8c-6080-11f1-ba76-6fd14a8cf585
md"""
# Spectral sensitivity with the Sh-matrix method

Fast, **exact** wavelength / refractive-index sensitivities for an axisymmetric
scatterer, built on the Sh-matrix moment-separation backend of
`TransitionMatrices.jl`.

The Sh-matrix method factors every EBCM surface integral as

```math
\int = \sum_{\text{terms}} \big[\,\text{coefficient}(k, m_r)\,\big] \times \big[\,\text{shape-only moment}\,\big],
```

where the **shape moments depend only on the geometry** — not on the wavelength
``\lambda`` or the refractive index ``m_r``. `prepare_sh` computes the moments
*once*; each `transition_matrix(prep, λ, mᵣ)` is then a cheap coefficient×moment
sum. Because that reconstruction is a plain differentiable function of
``(\lambda, m_r)`` with the moments riding through as constants, `ForwardDiff`
straight through it gives exact ``\partial T/\partial\lambda`` and
``\partial T/\partial m_r`` — at a fraction of the cost of differentiating the
full from-scratch assembly.
"""

# ╔═╡ 162c900e-6080-11f1-8b5d-bfd5a06ead20
md"""
## 1. Define the scatterer and prepare the moments

A prolate spheroid (semi-axes `a = 2`, `c = 1`) with a fixed small absorption
`mᵢ = 0.02`. `prepare_sh` runs the geometry quadrature **once** — every spectral
point below reuses it.
"""

# ╔═╡ 162c9018-6080-11f1-b436-23097d1541c2
begin
    nmax, Ng = 10, 200
    mᵢ = 0.02
    spheroid = TransitionMatrices.Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + mᵢ * im)
    prep = prepare_sh(spheroid, nmax, Ng)
end

# ╔═╡ 162c9022-6080-11f1-bb17-05963c0ec94f
md"""
## 2. Spectrum and sensitivities

Sweep the wavelength across a band; at each point reconstruct the T-matrix and
read off the scattering / extinction cross sections. The sensitivities
``\partial C_\text{sca}/\partial\lambda`` and ``\partial C_\text{sca}/\partial m_r``
come from `ForwardDiff` applied to the **reconstruction only** — the moments in
`prep` are constants.
"""

# ╔═╡ 162c9040-6080-11f1-b8d1-031603b7cffa
begin
    mᵣ = 1.5
    λs = collect(range(2π / 1.6, 2π / 0.6; length = 60))

    Csca(λ, m) = calc_Csca(transition_matrix(prep, λ, complex(m, mᵢ)), λ)
    Cext(λ, m) = calc_Cext(transition_matrix(prep, λ, complex(m, mᵢ)), λ)

    csca = [Csca(λ, mᵣ) for λ in λs]
    cext = [Cext(λ, mᵣ) for λ in λs]

    ∂Csca_∂λ = [ForwardDiff.derivative(l -> Csca(l, mᵣ), λ) for λ in λs]
    ∂Csca_∂mᵣ = [ForwardDiff.derivative(m -> Csca(λ, m), mᵣ) for λ in λs]
    nothing
end

# ╔═╡ 162c904a-6080-11f1-8f00-37f3c3cde468
md"""
## 3. The spectra and their sensitivities
"""

# ╔═╡ 162c9054-6080-11f1-99a7-5b1abfb33ebb
let
    p1 = plot(λs, [csca cext]; label = ["Csca" "Cext"], lw = 2,
        ylabel = "cross section", legend = :topright)
    p2 = plot(λs, ∂Csca_∂λ; label = "∂Csca/∂λ", lw = 2, color = 3, ylabel = "∂/∂λ")
    p3 = plot(λs, ∂Csca_∂mᵣ; label = "∂Csca/∂mᵣ", lw = 2, color = 4,
        xlabel = "wavelength λ", ylabel = "∂/∂mᵣ")
    plot(p1, p2, p3; layout = (3, 1), size = (720, 660),
        title = ["spheroid a=2, c=1, m=1.5+0.02im" "" ""])
end

# ╔═╡ 162c905e-6080-11f1-92ea-571ed46028d5
md"""
## 4. Why it is cheap — and how much

Differentiating the **reconstruction** reuses the single `prepare_sh` for the
whole spectrum; differentiating the **classic from-scratch assembly** re-runs the
full quadrature and Bessel recursions (in dual numbers) at every wavelength.
"""

# ╔═╡ 162c906a-6080-11f1-a118-e37de9a387e2
begin
    function Csca_classic(λ, m)
        T = promote_type(typeof(λ), typeof(m))
        s = TransitionMatrices.Spheroid{T, Complex{T}}(T(2.0), T(1.0), Complex{T}(m, mᵢ))
        calc_Csca(transition_matrix(s, λ, EBCM(nmax, Ng)), λ)
    end

    f_sh(λ) = ForwardDiff.derivative(l -> Csca(l, mᵣ), λ)
    f_classic(λ) = ForwardDiff.derivative(l -> Csca_classic(l, mᵣ), λ)
    f_sh(λs[1]);
    f_classic(λs[1])   # warm up

    t_prepare = @belapsed prepare_sh($spheroid, $nmax, $Ng) samples=1 evals=1
    t_sh = @belapsed [f_sh(λ) for λ in $λs] samples=1 evals=1
    t_classic = @belapsed [f_classic(λ) for λ in $λs] samples=1 evals=1

    (; n_points = length(λs),
        sh_ms = round((t_prepare + t_sh) * 1e3; digits = 1),
        classic_ms = round(t_classic * 1e3; digits = 1),
        speedup = round(t_classic / (t_prepare + t_sh); digits = 1))
end

# ╔═╡ 162c9072-6080-11f1-983e-97897cdba9f0
md"""
## Takeaway

`prepare_sh` once, then `ForwardDiff` the reconstruction, gives exact
``\partial C/\partial\lambda`` and ``\partial C/\partial m_r`` (matching the
classic result to machine precision) at a large speedup for spectra — the win
grows with the number of wavelengths because the geometry quadrature is
amortized.

Ideal for **material/wavelength sensitivity** (retrievals, dispersion fitting).
Shape-parameter derivatives act on the moments themselves and are better served
by the dedicated analytical EBCM linearization backend.
"""

# ╔═╡ Cell order:
# ╟─162c8f8c-6080-11f1-ba76-6fd14a8cf585
# ╠═162c8ff8-6080-11f1-af6b-f5508490c007
# ╟─162c900e-6080-11f1-8b5d-bfd5a06ead20
# ╠═162c9018-6080-11f1-b436-23097d1541c2
# ╟─162c9022-6080-11f1-bb17-05963c0ec94f
# ╠═162c9040-6080-11f1-b8d1-031603b7cffa
# ╟─162c904a-6080-11f1-8f00-37f3c3cde468
# ╠═162c9054-6080-11f1-99a7-5b1abfb33ebb
# ╟─162c905e-6080-11f1-92ea-571ed46028d5
# ╠═162c906a-6080-11f1-a118-e37de9a387e2
# ╟─162c9072-6080-11f1-983e-97897cdba9f0
