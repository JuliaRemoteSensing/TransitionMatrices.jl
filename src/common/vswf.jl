# ── Vector spherical wave functions (VSWFs) for near-field reconstruction ─────
#
# The incident, scattered and internal fields of a T-matrix problem are expanded
# in vector spherical wave functions. We use the convention of Mishchenko, Travis
# & Lacis (2002), Eqs. (5.16)–(5.17), which is the convention the rest of this
# package (`amplitude_matrix`, the T-matrix itself) already follows:
#
#   𝐌ₘₙ(k𝐫) = γₘₙ e^{imφ} zₙ(kr) [ i πₘₙ(ϑ) 𝛝̂ − τₘₙ(ϑ) 𝛗̂ ]
#   𝐍ₘₙ(k𝐫) = γₘₙ e^{imφ} { n(n+1)/(kr) zₙ(kr) Pₘₙ(ϑ) 𝐫̂
#                          + [kr·zₙ(kr)]′/(kr) [ τₘₙ(ϑ) 𝛝̂ + i πₘₙ(ϑ) 𝛗̂ ] }
#
# with
#   γₘₙ      = (−1)^m √[(2n+1) / (4π n(n+1))],
#   Pₘₙ(ϑ)  = d^n_{0m}(ϑ)         (Wigner d, = normalized associated Legendre),
#   πₘₙ(ϑ)  = m d^n_{0m}(ϑ)/sinϑ,
#   τₘₙ(ϑ)  = d d^n_{0m}(ϑ)/dϑ,
# and zₙ the spherical Bessel function jₙ (`kind = :regular`, used for the regular
# functions RgM/RgN that expand the incident and internal fields) or the spherical
# Hankel function hₙ⁽¹⁾ (`kind = :outgoing`, the radiating M/N that expand the
# scattered field, valid outside the circumscribing sphere).
#
# The radial factors are obtained from the package's Riccati–Bessel routines
# (`ricattibesselj`, `ricattibessely`):
#   jₙ(x) = ψₙ(x)/x,  yₙ(x) = χₙ(x)/x,  hₙ⁽¹⁾(x) = jₙ(x) + i yₙ(x),
#   [x·zₙ(x)]′/x = ψ′ₙ(x)/x      (regular)   resp.   (ψ′ₙ + i χ′ₙ)/x   (outgoing).
#
# 𝐌 and 𝐍 share the same γₘₙ, so they satisfy ∇×𝐌ₘₙ = k 𝐍ₘₙ and ∇×𝐍ₘₙ = k 𝐌ₘₙ
# regardless of the absolute value of γₘₙ; this curl identity is what the unit
# tests check. The absolute normalization γₘₙ is fixed to the Mishchenko value so
# that, combined with the plane-wave incident coefficients, the reconstructed
# field reproduces `amplitude_matrix` in the far zone (validated against Mie).

@inline function _sph_unit_vectors(ϑ::T, φ::T) where {T}
    sϑ, cϑ = sincos(ϑ)
    sφ, cφ = sincos(φ)
    r̂ = SVector(sϑ * cφ, sϑ * sφ, cϑ)
    ϑ̂ = SVector(cϑ * cφ, cϑ * sφ, -sϑ)
    φ̂ = SVector(-sφ, cφ, zero(T))
    return r̂, ϑ̂, φ̂
end

# Radial factors zₙ(x) and [x·zₙ(x)]′/x for n = 1:nₘₐₓ, plus zₙ(x)/x (the factor
# multiplying the 𝐫̂ component of 𝐍).
function _vswf_radial(kind::Symbol, nₘₐₓ::Integer, x::Number)
    if kind === :regular
        ψ, ψ′ = ricattibesselj(nₘₐₓ, estimate_ricattibesselj_extra_terms(nₘₐₓ, x), x)
        zₙ = ψ ./ x
        xzₙ′ = ψ′ ./ x
    elseif kind === :outgoing
        ψ, ψ′ = ricattibesselj(nₘₐₓ, estimate_ricattibesselj_extra_terms(nₘₐₓ, x), x)
        χ, χ′ = ricattibessely(nₘₐₓ, x)
        zₙ = (ψ .+ im .* χ) ./ x
        xzₙ′ = (ψ′ .+ im .* χ′) ./ x
    else
        throw(ArgumentError("`kind` must be :regular or :outgoing, got $kind"))
    end
    return zₙ, xzₙ′
end

@doc raw"""
```
vswf(kind::Symbol, m::Integer, n::Integer, k, r, ϑ, φ)
```

Evaluate the vector spherical wave functions ``\mathbf{M}_{mn}`` and
``\mathbf{N}_{mn}`` of degree `m` and order `n` at the point `(r, ϑ, φ)` for size
parameter `k` (so the radial argument is `x = k·r`), returning the pair `(𝐌, 𝐍)`
as `SVector{3}` of **Cartesian** field components.

`kind` selects the radial dependence:

- `:regular`  — spherical Bessel ``j_n`` (the regular `RgM`/`RgN`, finite at the
  origin; used for incident and internal fields).
- `:outgoing` — spherical Hankel ``h_n^{(1)}`` (the radiating `M`/`N` used for the
  scattered field; valid outside the circumscribing sphere).

The convention is Mishchenko, Travis & Lacis (2002), Eqs. (5.16)–(5.17),
consistent with [`amplitude_matrix`](@ref). `𝐌` and `𝐍` satisfy
``\nabla\times\mathbf{M}_{mn} = k\,\mathbf{N}_{mn}`` and
``\nabla\times\mathbf{N}_{mn} = k\,\mathbf{M}_{mn}``.
"""
function vswf(kind::Symbol, m::Integer, n::Integer, k::Real, r::Real, ϑ::Real, φ::Real)
    abs(m) ≤ n || throw(ArgumentError("require |m| ≤ n, got m=$m, n=$n"))
    T = float(promote_type(typeof(k), typeof(r), typeof(ϑ), typeof(φ)))
    ϑ, φ = T(ϑ), T(φ)
    x = T(k) * T(r)

    zₙ, xzₙ′ = _vswf_radial(kind, n, x)
    zn = zₙ[n]
    dzn = xzₙ′[n]
    zn_over_x = zn / x

    # Angular factors Pₘₙ = d^n_{0m}, τₘₙ = d(d^n_{0m})/dϑ, πₘₙ = m d^n_{0m}/sinϑ.
    dvec, τvec = wigner_d_recursion(T, 0, m, n, ϑ; deriv = true)
    P = dvec[n]
    τ = τvec[n]
    π_ = pi_func(T, m, n, ϑ; d = P)

    γ = (-1)^(m & 1) * sqrt(T(2n + 1) / (4 * T(π) * n * (n + 1)))
    eⁱᵐᵠ = cis(m * φ)
    r̂, ϑ̂, φ̂ = _sph_unit_vectors(ϑ, φ)

    pref = γ * eⁱᵐᵠ
    𝐌 = pref * zn * (im * π_ .* ϑ̂ .- τ .* φ̂)
    𝐍 = pref *
        (n * (n + 1) * zn_over_x * P .* r̂ .+ dzn .* (τ .* ϑ̂ .+ im * π_ .* φ̂))

    return 𝐌, 𝐍
end

"""
```
vswf_cartesian(kind::Symbol, m, n, k, xyz::AbstractVector)
```

Convenience wrapper of [`vswf`](@ref) taking a Cartesian position `xyz = (x, y, z)`
and returning the `(𝐌, 𝐍)` Cartesian field vectors. The polar axis is `+z`.
"""
function vswf_cartesian(kind::Symbol, m::Integer, n::Integer, k::Real,
        xyz::AbstractVector)
    x, y, z = xyz[1], xyz[2], xyz[3]
    r = sqrt(x^2 + y^2 + z^2)
    ϑ = acos(clamp(z / r, -one(r), one(r)))
    φ = atan(y, x)
    return vswf(kind, m, n, k, r, ϑ, φ)
end

@testitem "VSWF satisfy the curl identities ∇×M = kN and ∇×N = kM" begin
    using TransitionMatrices: vswf_cartesian

    # Finite-difference Cartesian curl of a vector field F(𝐫).
    function fd_curl(F, p, h)
        ex = [h, 0.0, 0.0]
        ey = [0.0, h, 0.0]
        ez = [0.0, 0.0, h]
        ∂x = (F(p + ex) - F(p - ex)) / 2h
        ∂y = (F(p + ey) - F(p - ey)) / 2h
        ∂z = (F(p + ez) - F(p - ez)) / 2h
        return [∂y[3] - ∂z[2], ∂z[1] - ∂x[3], ∂x[2] - ∂y[1]]
    end

    k = 1.3
    h = 1e-5
    # A handful of generic off-axis points (avoid the poles / origin).
    points = ([0.7, 0.5, 0.9], [-0.6, 0.8, 0.4], [1.1, -0.3, -0.7])

    @testset "kind = $kind" for kind in (:regular, :outgoing)
        for p in points
            for n in 1:4, m in (-n):n
                M(q) = vswf_cartesian(kind, m, n, k, q)[1]
                N(q) = vswf_cartesian(kind, m, n, k, q)[2]
                curlM = fd_curl(M, p, h)
                curlN = fd_curl(N, p, h)
                kN = k * vswf_cartesian(kind, m, n, k, p)[2]
                kM = k * vswf_cartesian(kind, m, n, k, p)[1]
                @test maximum(abs, curlM - kN) / maximum(abs, kN) < 1e-6
                @test maximum(abs, curlN - kM) / maximum(abs, kM) < 1e-6
            end
        end
    end
end
