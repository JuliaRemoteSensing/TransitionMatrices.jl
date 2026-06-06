### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 7a1f0001-0000-4000-8000-000000000001
begin
    import Pkg
    Pkg.activate(@__DIR__)                 # the examples/ environment
    Pkg.develop(; path = dirname(@__DIR__)) # point TransitionMatrices at this repo
    Pkg.instantiate()
    using TransitionMatrices, Plots
end

# ╔═╡ 7a1f0002-0000-4000-8000-000000000002
md"""
# Near-field maps from a T-matrix

`TransitionMatrices.jl` reconstructs the electromagnetic field *around* a particle
— not just the far-field cross sections — from any computed T-matrix. Given an
incident plane wave, the scattered field is expanded in radiating vector spherical
wave functions (VSWFs) with coefficients ``(p, q) = \mathbf{T}\,(a, b)``, where
``(a, b)`` are the plane-wave expansion coefficients. The total external field is

```math
\mathbf{E}_\text{tot}(\mathbf r) = \mathbf{E}_\text{inc}(\mathbf r) + \mathbf{E}_\text{sca}(\mathbf r),
\qquad
\mathbf{E}_\text{sca}(\mathbf r) = \sum_{n,m} p_{mn}\mathbf{M}_{mn}(k\mathbf r) + q_{mn}\mathbf{N}_{mn}(k\mathbf r).
```

This works for **any** T-matrix (Mie, EBCM, IITM, Sh-matrix). Here we map the field
enhancement ``|\mathbf{E}_\text{tot}|/|\mathbf{E}_\text{inc}|`` around a dielectric
sphere.

!!! note "Region of validity"
    The radiating expansion converges only **outside the sphere circumscribing the
    particle** (the Rayleigh hypothesis). For a sphere that boundary *is* the
    surface, so the map is valid right down to it; for a non-spherical particle,
    evaluate only outside the circumscribing sphere.
"""

# ╔═╡ 7a1f0003-0000-4000-8000-000000000003
md"""
## 1. A scatterer and its incidence

A dielectric sphere of size parameter ``x = k a = 2.5`` and refractive index
``m_r = 1.5``, illuminated by an ``x``-polarized plane wave propagating along
``+z``. With ``\lambda = 2\pi`` the size parameter equals the radius, so ``a = 2.5``.

The incidence direction is ``(\vartheta_\text{inc}, \varphi_\text{inc}) = (0, 0)``;
the polarization ``(E_\vartheta, E_\varphi) = (1, 0)`` is ``\hat{\mathbf x}`` because
``\hat{\boldsymbol\vartheta}(0,0) = \hat{\mathbf x}``.
"""

# ╔═╡ 7a1f0004-0000-4000-8000-000000000004
begin
    λ = 2π
    k = 2π / λ
    a = 2.5
    mᵣ = 1.5 + 0.0im
    x = k * a
    N = ceil(Int, x + 4 * cbrt(x) + 2)
    𝐓 = TransitionMatrices.MieTransitionMatrix{ComplexF64, N}(x, mᵣ)

    # Incident plane wave: +z propagation, x-polarized.
    ϑ_inc, φ_inc = 0.0, 0.0
    Eθ, Eφ = 1.0 + 0im, 0.0im

    # The scattered coefficients (p, q) = 𝐓 (a, b) depend only on the incidence,
    # so compute them once and reuse them at every field point.
    p, q = scattering_coefficients(𝐓, ϑ_inc, φ_inc, Eθ, Eφ)
    nothing
end

# ╔═╡ 7a1f0005-0000-4000-8000-000000000005
md"""
## 2. Sample the total field on a grid

We evaluate ``|\mathbf{E}_\text{tot}|`` in the ``x``–``z`` plane (the plane of
incidence and polarization). Points inside the sphere are masked — the external
expansion is not meant to represent the interior field (that is the *internal*
field, a separate reconstruction).
"""

# ╔═╡ 7a1f0006-0000-4000-8000-000000000006
begin
    xs = range(-3a, 3a; length = 161)
    zs = range(-3a, 3a; length = 161)
    enhancement = fill(NaN, length(zs), length(xs))
    for (i, zz) in enumerate(zs), (j, xx) in enumerate(xs)
        r = hypot(xx, zz)
        r ≤ a && continue                      # inside the sphere: masked
        pos = [xx, 0.0, zz]
        E = incident_field(λ, ϑ_inc, φ_inc, Eθ, Eφ, pos) +
            scattered_field(p, q, λ, pos)
        enhancement[i, j] = sqrt(abs2(E[1]) + abs2(E[2]) + abs2(E[3]))
    end
    nothing
end

# ╔═╡ 7a1f0007-0000-4000-8000-000000000007
md"""
## 3. The field-enhancement map

The incident wave travels upward (``+z``). The map shows the standing-wave
interference of the incident and scattered fields, the forward-scattering
concentration on the far (``+z``) side, and the shadow on the near side.
"""

# ╔═╡ 7a1f0008-0000-4000-8000-000000000008
let
    hm = heatmap(xs, zs, enhancement;
        aspect_ratio = 1, c = :inferno, clims = (0, maximum(filter(isfinite, enhancement))),
        xlabel = "x", ylabel = "z", colorbar_title = "|E| / |E₀|",
        title = "dielectric sphere, x = $(round(x; digits = 2)), mᵣ = 1.5",
        size = (640, 560))
    # outline the sphere
    θ = range(0, 2π; length = 200)
    plot!(hm, a .* cos.(θ), a .* sin.(θ); label = "", lw = 1.5, lc = :cyan)
end

# ╔═╡ 7a1f0009-0000-4000-8000-000000000009
md"""
## 4. A line cut along the optical axis

Along ``x = 0`` the field shows the shadow/illumination asymmetry and the
forward enhancement directly.
"""

# ╔═╡ 7a1f000a-0000-4000-8000-00000000000a
let
    zline = range(-3a, 3a; length = 400)
    amp = map(zline) do zz
        abs(zz) ≤ a && return NaN
        pos = [0.0, 0.0, zz]
        E = incident_field(λ, ϑ_inc, φ_inc, Eθ, Eφ, pos) + scattered_field(p, q, λ, pos)
        sqrt(abs2(E[1]) + abs2(E[2]) + abs2(E[3]))
    end
    plot(zline, amp; lw = 2, xlabel = "z (optical axis)", ylabel = "|E| / |E₀|",
        label = "", title = "axial cut (x = y = 0)", size = (720, 320))
    vspan!([-a, a]; alpha = 0.15, c = :gray, label = "sphere")
end

# ╔═╡ 7a1f000b-0000-4000-8000-00000000000b
md"""
## Takeaway

`scattering_coefficients` once, then `scattered_field` / `total_field` at any point
gives the full external field for any T-matrix — sphere, spheroid, cylinder,
Chebyshev, or a general particle. Reusing the precomputed `(p, q)` keeps a dense
field grid cheap.

For the field **inside** the particle, see the internal-field reconstruction
(`TODO`, Tier 2).
"""

# ╔═╡ Cell order:
# ╟─7a1f0002-0000-4000-8000-000000000002
# ╠═7a1f0001-0000-4000-8000-000000000001
# ╟─7a1f0003-0000-4000-8000-000000000003
# ╠═7a1f0004-0000-4000-8000-000000000004
# ╟─7a1f0005-0000-4000-8000-000000000005
# ╠═7a1f0006-0000-4000-8000-000000000006
# ╟─7a1f0007-0000-4000-8000-000000000007
# ╠═7a1f0008-0000-4000-8000-000000000008
# ╟─7a1f0009-0000-4000-8000-000000000009
# ╠═7a1f000a-0000-4000-8000-00000000000a
# ╟─7a1f000b-0000-4000-8000-00000000000b
