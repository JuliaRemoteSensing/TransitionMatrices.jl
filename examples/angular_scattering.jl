### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 226c9e5c-6082-11f1-9059-f9e778e9b2f8
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.develop(; path = dirname(@__DIR__))
    Pkg.instantiate()
    using TransitionMatrices, Plots
end

# ╔═╡ 226c9e04-6082-11f1-9773-191ee2516176
md"""
# Angular scattering & polarization

From a T-matrix you can read off the full far field: the **amplitude matrix**
``\mathbf{S}`` at a fixed geometry, its **phase matrix** ``\mathbf{Z}``, and the
orientation-averaged **scattering (Mueller) matrix** ``\mathbf{F}(\theta)``.

This notebook plots the phase function ``F_{11}(\theta)``, the degree of linear
polarization ``-F_{12}/F_{11}`` (for unpolarized incidence) and the ratios
``F_{22}/F_{11}``, ``F_{44}/F_{11}`` for a spheroid.
"""

# ╔═╡ 226c9e70-6082-11f1-a65b-3b159155fbf8
md"""
## Build a T-matrix
"""

# ╔═╡ 226c9e7a-6082-11f1-9dc0-75a4e1f6bbbb
begin
    λ = 2π
    spheroid = TransitionMatrices.Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)
    T = calc_T(spheroid, λ)
end

# ╔═╡ 226c9e84-6082-11f1-bc9f-71d54665c33b
md"""
## Orientation-averaged scattering matrix vs angle

`scattering_matrix(T, λ, θs)` returns an `Nθ × 6` matrix whose columns are
`[F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄]` (`θs` in degrees).
"""

# ╔═╡ 226c9e8e-6082-11f1-bdcf-7b0495b58d1e
begin
    θs = collect(0.0:1.0:180.0)
    F = scattering_matrix(T, λ, θs)
    F11 = F[:, 1]
    pol = -F[:, 2] ./ F[:, 1]      # degree of linear polarization
    F22_11 = F[:, 3] ./ F[:, 1]
    F44_11 = F[:, 6] ./ F[:, 1]
    nothing
end

# ╔═╡ 226c9eac-6082-11f1-8892-375f3d9f7416
md"""
## Plots
"""

# ╔═╡ 226c9eb6-6082-11f1-9cb8-e38b0c25c94a
let
    p1 = plot(θs, F11; yscale = :log10, lw = 2, label = "F₁₁",
        ylabel = "phase function", legend = :topright)
    p2 = plot(θs, pol; lw = 2, color = 2, label = "-F₁₂/F₁₁",
        ylabel = "lin. polarization")
    p3 = plot(θs, [F22_11 F44_11]; lw = 2, label = ["F₂₂/F₁₁" "F₄₄/F₁₁"],
        xlabel = "scattering angle θ (deg)", ylabel = "ratio")
    plot(p1, p2, p3; layout = (3, 1), size = (720, 660),
        title = ["spheroid a=2, c=1, m=1.5+0.02im" "" ""])
end

# ╔═╡ 226c9eca-6082-11f1-a933-f9b57543d7a2
md"""
## Amplitude and phase matrix at a fixed geometry

For a specific incidence/scattering direction, `amplitude_matrix` gives the 2×2
complex ``\mathbf{S}`` and `phase_matrix` the 4×4 Mueller ``\mathbf{Z}``.
"""

# ╔═╡ 226c9ede-6082-11f1-82d2-43db7f93d579
begin
    S = amplitude_matrix(T, 0.0, 0.0, π / 4, 0.0; λ = λ)   # ϑᵢ,φᵢ → ϑₛ,φₛ
    Z = phase_matrix(S)
    (; S, Z)
end

# ╔═╡ Cell order:
# ╟─226c9e04-6082-11f1-9773-191ee2516176
# ╠═226c9e5c-6082-11f1-9059-f9e778e9b2f8
# ╟─226c9e70-6082-11f1-a65b-3b159155fbf8
# ╠═226c9e7a-6082-11f1-9dc0-75a4e1f6bbbb
# ╟─226c9e84-6082-11f1-bc9f-71d54665c33b
# ╠═226c9e8e-6082-11f1-bdcf-7b0495b58d1e
# ╟─226c9eac-6082-11f1-8892-375f3d9f7416
# ╠═226c9eb6-6082-11f1-9cb8-e38b0c25c94a
# ╟─226c9eca-6082-11f1-a933-f9b57543d7a2
# ╠═226c9ede-6082-11f1-82d2-43db7f93d579
