# ── Internal field (Tier 2): general axisymmetric particle (EBCM) ─────────────
#
# Extends `internal_coefficients` / `internal_field` (defined for the homogeneous
# sphere in `common/near_field.jl`) to a general axisymmetric particle. The
# interior coefficients come from the EBCM matrices rather than an analytic
# formula. With the package's 𝐐 = 𝐏 + i𝐔 (so 𝐓 = -𝐏 𝐐⁻¹) the Waterman relation
# gives the internal coefficients as 𝐐⁻¹ 𝐚 per azimuthal block; an overall factor
# ½ — universal across size and refractive index, fixed by matching the analytic
# Mie internal field on a sphere to machine precision — converts to this package's
# VSWF normalization:
#
#   𝐜⁽ᵐ⁾ = ½ 𝐐⁽ᵐ⁾⁻¹ 𝐚⁽ᵐ⁾      (degree +m).
#
# The −m block follows the same ±m symmetry the T-matrix uses (see
# `AxisymmetricTransitionMatrix`): the M↔N cross blocks flip sign, so
# 𝐜⁽⁻ᵐ⁾ = 𝐃 (½ 𝐐⁽ᵐ⁾⁻¹) 𝐃 𝐚⁽⁻ᵐ⁾ with 𝐃 = diag(+𝐈_M, -𝐈_N). Each block stacks the
# M (p=1) coefficients above the N (p=2) coefficients over n = max(1,m):nmax.
#
# VALIDITY: the single-origin regular-VSWF interior expansion converges only
# WITHIN THE INSCRIBED SPHERE of the particle (radius `rmin(shape)`); between it
# and the surface the series may diverge. For a sphere the inscribed sphere is the
# whole interior. The convention is validated on the degenerate sphere (a spheroid
# with a = c) against the analytic Mie internal field; for genuinely non-spherical
# shapes there is no self-contained surface check (tangential boundary continuity
# would need both the interior and exterior expansions to converge at the surface,
# which the Rayleigh hypothesis forbids).

@doc raw"""
```
internal_coefficients(shape::AbstractAxisymmetricShape, λ, nmax, Ng, ϑ_inc, φ_inc, Eθ, Eφ) -> (c, d)
```

Internal-field expansion coefficients `(cₘₙ, dₘₙ)` for a general axisymmetric
particle, obtained from the EBCM matrices: per azimuthal block ``\mathbf{c} =
\tfrac{1}{2}\mathbf{Q}^{-1}\mathbf{a}`` with ``\mathbf{Q} = \mathbf{P} +
\mathrm{i}\mathbf{U}``, under an incident plane wave (`(ϑ_inc, φ_inc)` propagation,
Jones polarization `(Eθ, Eφ)`). Returns `OffsetArray`s indexed `[m, n]`; reuse them
across field points with [`internal_field`](@ref).

The reconstruction is valid only inside the particle's **inscribed sphere** (radius
`rmin(shape)`).
"""
function internal_coefficients(shape::AbstractAxisymmetricShape{T, CT}, λ, nmax::Integer,
        Ng::Integer, ϑ_inc, φ_inc, Eθ, Eφ) where {T, CT}
    Tr = real(CT)
    a, b = _plane_wave_coefficients(nmax, Tr(ϑ_inc), Tr(φ_inc), CT(Eθ), CT(Eφ))
    c = OffsetArray(zeros(CT, 2nmax + 1, nmax), (-nmax):nmax, 1:nmax)
    d = OffsetArray(zeros(CT, 2nmax + 1, nmax), (-nmax):nmax, 1:nmax)
    for mₐ in 0:nmax
        𝐏, 𝐔 = mₐ == 0 ? ebcm_matrices_m₀(shape, λ, nmax, Ng) :
                ebcm_matrices_m(mₐ, shape, λ, nmax, Ng)
        𝐐 = 𝐏 .+ im .* 𝐔
        ns = max(1, mₐ):nmax
        nn = length(ns)
        avec = CT[[a[mₐ, n] for n in ns]; [b[mₐ, n] for n in ns]]
        cvec = (𝐐 \ avec) ./ 2
        for (i, n) in enumerate(ns)
            c[mₐ, n] = cvec[i]
            d[mₐ, n] = cvec[nn + i]
        end
        if mₐ > 0
            𝐃 = CT[ones(CT, nn); -ones(CT, nn)]
            amvec = CT[[a[-mₐ, n] for n in ns]; [b[-mₐ, n] for n in ns]]
            cmvec = (𝐃 .* (𝐐 \ (𝐃 .* amvec))) ./ 2
            for (i, n) in enumerate(ns)
                c[-mₐ, n] = cmvec[i]
                d[-mₐ, n] = cmvec[nn + i]
            end
        end
    end
    return c, d
end

@doc raw"""
```
internal_field(shape::AbstractAxisymmetricShape, λ, nmax, Ng, ϑ_inc, φ_inc, Eθ, Eφ, r⃗) -> SVector{3}
```

Internal electric field inside a general axisymmetric particle at the Cartesian
point `r⃗`, via the EBCM internal coefficients (see [`internal_coefficients`](@ref)).
The relative refractive index is taken from `shape`. `r⃗` must lie within the
particle's **inscribed sphere** (radius `rmin(shape)`), where the interior
expansion converges.
"""
function internal_field(shape::AbstractAxisymmetricShape, λ, nmax::Integer, Ng::Integer,
        ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
    c, d = internal_coefficients(shape, λ, nmax, Ng, ϑ_inc, φ_inc, Eθ, Eφ)
    return internal_field(c, d, shape.m, λ, r⃗)
end

@testitem "EBCM internal field matches Mie on a sphere (degenerate spheroid)" begin
    using TransitionMatrices: Spheroid, internal_field

    # A sphere is a spheroid with a = c. Its EBCM internal field must reproduce the
    # analytic Mie internal field — this locks the EBCM internal-coefficient
    # convention (the ½ factor and the ±m symmetry) to machine precision.
    λ = 2π
    k0 = 2π / λ
    @testset "R=$R, mᵣ=$mᵣ" for (R, mᵣ) in ((2.0, 1.5 + 0.05im), (1.5, 1.33 + 0.0im))
        x = k0 * R
        nmax = ceil(Int, x + 4cbrt(x) + 2)
        Ng = 200
        sph = Spheroid{Float64, ComplexF64}(R, R, mᵣ)
        worst = 0.0
        for ϑ in (0.4, 1.1, 2.0), φ in (0.0, 1.3, 4.0), rr in (0.3R, 0.6R, 0.85R)
            sϑ, cϑ, sφ, cφ = sin(ϑ), cos(ϑ), sin(φ), cos(φ)
            pos = [rr * sϑ * cφ, rr * sϑ * sφ, rr * cϑ]
            Eebcm = internal_field(sph, λ, nmax, Ng, 0.0, 0.0, 1.0 + 0im, 0.0im, pos)
            Emie = internal_field(x, mᵣ, λ, 0.0, 0.0, 1.0 + 0im, 0.0im, pos; nmax)
            worst = max(worst, maximum(abs, Eebcm - Emie) / maximum(abs, Emie))
        end
        @test worst < 1e-8
    end
end

@testitem "EBCM internal field reconstructs inside a spheroid (runs, finite)" begin
    using TransitionMatrices: Spheroid, internal_coefficients, internal_field

    # Non-spherical: no self-contained reference, but the coefficients/field must be
    # finite and the two-step and one-shot forms must agree, inside the inscribed
    # sphere (radius = min semi-axis).
    λ = 2π
    sph = Spheroid{Float64, ComplexF64}(1.0, 2.0, 1.4 + 0.02im)  # prolate, rmin = 1
    nmax, Ng = 10, 200
    c, d = internal_coefficients(sph, λ, nmax, Ng, 0.3, 0.5, 1.0 + 0im, 0.0im)
    pos = [0.4, 0.2, 0.5]                       # |pos| ≈ 0.67 < inscribed radius 1
    E1 = internal_field(c, d, sph.m, λ, pos)
    E2 = internal_field(sph, λ, nmax, Ng, 0.3, 0.5, 1.0 + 0im, 0.0im, pos)
    @test all(isfinite, abs.(E1))
    @test maximum(abs, E1 - E2) < 1e-12
end
