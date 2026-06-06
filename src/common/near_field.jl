# ── Near-field reconstruction (Tier 1: external field) ────────────────────────
#
# Given a T-matrix `𝐓` and an incident plane wave, reconstruct the electric field
# at points OUTSIDE the circumscribing sphere of the particle: the incident field
# (analytic plane wave), the scattered field, and their sum (the total field).
#
# The incident plane wave propagates along the direction `n̂ = (ϑ_inc, φ_inc)` and
# carries a polarization given as a Jones vector `(Eθ, Eφ)` in the spherical basis
# `(𝛝̂, 𝛗̂)` evaluated at `n̂`:
#
#   𝐄_inc(𝐫) = (Eθ 𝛝̂ + Eφ 𝛗̂) exp(i k n̂·𝐫).
#
# It is expanded in regular VSWFs (`vswf(:regular, …)`)
#   𝐄_inc(𝐫) = Σ_{n,m} [ aₘₙ RgMₘₙ(k𝐫) + bₘₙ RgNₘₙ(k𝐫) ],
# with coefficients (locked numerically against the analytic plane wave to machine
# precision, Mishchenko convention; γₘₙ as in `vswf`):
#   aₘₙ = 4π iⁿ   γₘₙ e^{-imφ_inc} ( −i πₘₙ(ϑ_inc) Eθ − τₘₙ(ϑ_inc) Eφ ),
#   bₘₙ = 4π iⁿ⁻¹ γₘₙ e^{-imφ_inc} (   τₘₙ(ϑ_inc) Eθ − i πₘₙ(ϑ_inc) Eφ ).
#
# The scattered field uses the radiating VSWFs and the scattered coefficients
# (p,q) = 𝐓 (a,b):
#   𝐄_sca(𝐫) = Σ_{n,m} [ pₘₙ Mₘₙ(k𝐫) + qₘₙ Nₘₙ(k𝐫) ],
#   pₘₙ = Σ_{n′m′} 𝐓₁₁ aₘ′ₙ′ + 𝐓₁₂ bₘ′ₙ′,   qₘₙ = Σ_{n′m′} 𝐓₂₁ aₘ′ₙ′ + 𝐓₂₂ bₘ′ₙ′.
#
# Validity: the radiating expansion converges only OUTSIDE the smallest sphere
# circumscribing the particle (the Rayleigh hypothesis). Inside that sphere the
# scattered series may diverge; callers are responsible for evaluating at r larger
# than the circumscribing radius.

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

@inline function _cart_to_sph(r⃗)
    x, y, z = r⃗[1], r⃗[2], r⃗[3]
    r = sqrt(x^2 + y^2 + z^2)
    ϑ = acos(clamp(z / r, -one(r), one(r)))
    φ = atan(y, x)
    return r, ϑ, φ
end

# Plane-wave expansion coefficients aₘₙ, bₘₙ (degree m, order n), as OffsetArrays
# indexed `[m, n]`, m ∈ -N:N, n ∈ 1:N.
function _plane_wave_coefficients(N::Integer, ϑ::T, φ::T, Eθ::Complex{T},
        Eφ::Complex{T}) where {T}
    a = OffsetArray(zeros(Complex{T}, 2N + 1, N), (-N):N, 1:N)
    b = OffsetArray(zeros(Complex{T}, 2N + 1, N), (-N):N, 1:N)
    for m in (-N):N
        dvec, τvec = wigner_d_recursion(T, 0, m, N, ϑ; deriv = true)
        for n in max(1, abs(m)):N
            P = dvec[n]
            τ = τvec[n]
            π_ = pi_func(T, m, n, ϑ; d = P)
            γ = (-1)^(m & 1) * sqrt(T(2n + 1) / (4 * T(π) * n * (n + 1)))
            pref = 4 * T(π) * cis(-m * φ) * γ
            a[m, n] = pref * im^(n & 3) * (-im * π_ * Eθ - τ * Eφ)
            b[m, n] = pref * im^((n - 1) & 3) * (τ * Eθ - im * π_ * Eφ)
        end
    end
    return a, b
end

@doc raw"""
```
scattering_coefficients(𝐓, ϑ_inc, φ_inc, Eθ, Eφ) -> (p, q)
```

Scattered-field expansion coefficients `(pₘₙ, qₘₙ)` for the T-matrix `𝐓` under an
incident plane wave propagating along `(ϑ_inc, φ_inc)` with Jones polarization
`(Eθ, Eφ)` in the spherical basis at the incidence direction. Returns two
`OffsetArray`s indexed `[m, n]` (m ∈ -N:N, n ∈ 1:N), where `(p, q) = 𝐓 (a, b)` and
`(a, b)` are the incident plane-wave coefficients. Reuse these across many field
points (see [`scattered_field`](@ref)).
"""
function scattering_coefficients(𝐓::AbstractTransitionMatrix{CT, N}, ϑ_inc, φ_inc,
        Eθ, Eφ) where {CT, N}
    T = real(CT)
    a, b = _plane_wave_coefficients(N, T(ϑ_inc), T(φ_inc), CT(Eθ), CT(Eφ))
    p = OffsetArray(zeros(CT, 2N + 1, N), (-N):N, 1:N)
    q = OffsetArray(zeros(CT, 2N + 1, N), (-N):N, 1:N)
    for n in 1:N, m in (-n):n
        s1 = zero(CT)
        s2 = zero(CT)
        for n′ in 1:N, m′ in (-n′):n′
            aₘ′ₙ′ = a[m′, n′]
            bₘ′ₙ′ = b[m′, n′]
            s1 += 𝐓[m, n, m′, n′, 1, 1] * aₘ′ₙ′ + 𝐓[m, n, m′, n′, 1, 2] * bₘ′ₙ′
            s2 += 𝐓[m, n, m′, n′, 2, 1] * aₘ′ₙ′ + 𝐓[m, n, m′, n′, 2, 2] * bₘ′ₙ′
        end
        p[m, n] = s1
        q[m, n] = s2
    end
    return p, q
end

# Sum Σ c₁ₘₙ Fₘₙ + c₂ₘₙ Gₘₙ at one point, where (F,G)=(RgM,RgN) for kind=:regular
# or (M,N) for kind=:outgoing. Radial factors computed once; angular per m.
function _field_from_coeffs(c₁, c₂, kind::Symbol, k::Real, r⃗, N::Integer)
    T = float(promote_type(typeof(k), eltype(r⃗)))
    r, ϑ, φ = _cart_to_sph(SVector{3, T}(r⃗[1], r⃗[2], r⃗[3]))
    x = T(k) * r
    zₙ, xzₙ′ = _vswf_radial(kind, N, x)
    r̂, ϑ̂, φ̂ = _sph_unit_vectors(ϑ, φ)

    CT = complex(T)
    E = zero(SVector{3, CT})
    for m in (-N):N
        dvec, τvec = wigner_d_recursion(T, 0, m, N, ϑ; deriv = true)
        eⁱᵐᵠ = cis(m * φ)
        for n in max(1, abs(m)):N
            P = dvec[n]
            τ = τvec[n]
            π_ = pi_func(T, m, n, ϑ; d = P)
            γ = (-1)^(m & 1) * sqrt(T(2n + 1) / (4 * T(π) * n * (n + 1)))
            pref = γ * eⁱᵐᵠ
            zn = zₙ[n]
            dzn = xzₙ′[n]
            znx = zn / x
            𝐅 = pref * zn * (im * π_ .* ϑ̂ .- τ .* φ̂)
            𝐆 = pref *
                (n * (n + 1) * znx * P .* r̂ .+ dzn .* (τ .* ϑ̂ .+ im * π_ .* φ̂))
            E += c₁[m, n] * 𝐅 + c₂[m, n] * 𝐆
        end
    end
    return E
end

@doc raw"""
```
incident_field(λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗) -> SVector{3}
```

The incident plane wave ``(E_\theta\hat{\boldsymbol\vartheta} + E_\varphi\hat{\boldsymbol\varphi})\exp(\mathrm{i}k\,\hat{\mathbf n}\cdot\mathbf r)``
evaluated at the Cartesian point `r⃗`, propagating along `n̂ = (ϑ_inc, φ_inc)` with
`k = 2π/λ`. Polarization `(Eθ, Eφ)` is given in the spherical basis at `n̂`.
"""
function incident_field(λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
    T = float(promote_type(typeof(λ), typeof(ϑ_inc), typeof(φ_inc), eltype(r⃗)))
    k = 2 * T(π) / T(λ)
    sϑ, cϑ = sincos(T(ϑ_inc))
    sφ, cφ = sincos(T(φ_inc))
    n̂ = SVector(sϑ * cφ, sϑ * sφ, cϑ)
    ϑ̂ = SVector(cϑ * cφ, cϑ * sφ, -sϑ)
    φ̂ = SVector(-sφ, cφ, zero(T))
    return (Complex{T}(Eθ) .* ϑ̂ .+ Complex{T}(Eφ) .* φ̂) .* cis(k * _dot3(n̂, r⃗))
end

@doc raw"""
```
scattered_field(p, q, λ, r⃗) -> SVector{3}
scattered_field(𝐓, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗) -> SVector{3}
```

Scattered electric field at the Cartesian point `r⃗`, reconstructed from the
radiating VSWFs. The first form reuses precomputed coefficients `(p, q)` from
[`scattering_coefficients`](@ref) (efficient for evaluating many points); the
second computes them on the fly for the T-matrix `𝐓` and incidence.

The point `r⃗` must lie OUTSIDE the sphere circumscribing the particle (the
radiating series diverges inside it).
"""
function scattered_field(p::OffsetArray, q::OffsetArray, λ, r⃗)
    N = size(p, 2)
    Tr = real(eltype(p))
    k = 2 * Tr(π) / λ
    return _field_from_coeffs(p, q, :outgoing, k, r⃗, N)
end

function scattered_field(𝐓::AbstractTransitionMatrix, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
    p, q = scattering_coefficients(𝐓, ϑ_inc, φ_inc, Eθ, Eφ)
    return scattered_field(p, q, λ, r⃗)
end

@doc raw"""
```
total_field(𝐓, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗) -> SVector{3}
```

Total external field `𝐄_inc(r⃗) + 𝐄_sca(r⃗)` at the Cartesian point `r⃗` for the
T-matrix `𝐓` under the incident plane wave (see [`incident_field`](@ref) and
[`scattered_field`](@ref)).
"""
function total_field(𝐓::AbstractTransitionMatrix, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
    return incident_field(λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗) +
           scattered_field(𝐓, λ, ϑ_inc, φ_inc, Eθ, Eφ, r⃗)
end

@testitem "Incident plane-wave expansion reconstructs the analytic plane wave" begin
    using TransitionMatrices: _plane_wave_coefficients, _field_from_coeffs, incident_field

    # Σ aₘₙ RgMₘₙ + bₘₙ RgNₘₙ must equal the analytic plane wave inside the
    # circumscribing sphere — this locks the incident coefficients and the regular
    # VSWF normalization simultaneously.
    λ = 2π
    k = 2π / λ
    N = 20
    pts = ([0.2, 0.1, 0.15], [-0.18, 0.25, -0.2], [0.05, -0.3, 0.22])
    for (Eθ, Eφ) in ((1.0 + 0im, 0.0im), (0.0im, 1.0 + 0im), (0.6 + 0.2im, 0.3 - 0.4im)),
        (ϑi, φi) in ((0.7, 0.4), (2.1, -1.3))

        a, b = _plane_wave_coefficients(N, ϑi, φi, ComplexF64(Eθ), ComplexF64(Eφ))
        for p in pts
            Erec = _field_from_coeffs(a, b, :regular, k, p, N)
            Epw = incident_field(λ, ϑi, φi, Eθ, Eφ, p)
            @test maximum(abs, Erec - Epw) < 1e-10
        end
    end
end

@testitem "Scattered field reproduces amplitude_matrix in the far zone" begin
    using TransitionMatrices: MieTransitionMatrix, scattering_coefficients,
                              scattered_field, amplitude_matrix

    dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    function bvecs(ϑ, φ)
        sϑ, cϑ = sincos(ϑ)
        sφ, cφ = sincos(φ)
        [sϑ * cφ, sϑ * sφ, cϑ], [cϑ * cφ, cϑ * sφ, -sϑ], [-sφ, cφ, 0.0]
    end

    λ = 2π
    k = 2π / λ
    x = 3.0
    N = ceil(Int, x + 4cbrt(x) + 2)
    𝐓 = MieTransitionMatrix{ComplexF64, N}(x, 1.5 + 0.01im)
    Eθ, Eφ = 0.5 + 0.3im, -0.2 + 0.4im
    ϑi, φi = 0.6, 0.5
    p, q = scattering_coefficients(𝐓, ϑi, φi, Eθ, Eφ)

    # The reconstructed far field converges to (e^{ikR}/R)·S·E₀ as O(1/R); checking
    # the ratio at two radii a decade apart confirms both the value and the rate,
    # which locks the outgoing-VSWF normalization and the (p,q)=𝐓(a,b) application.
    for (ϑs, φs) in ((0.4, 0.3), (1.5, 2.0), (2.6, -0.7))
        n̂s, θ̂s, φ̂s = bvecs(ϑs, φs)
        S = amplitude_matrix(𝐓, ϑi, φi, ϑs, φs; λ)
        predθ(R) = cis(k * R) / R * (S[1, 1] * Eθ + S[1, 2] * Eφ)
        predφ(R) = cis(k * R) / R * (S[2, 1] * Eθ + S[2, 2] * Eφ)
        e = Float64[]
        for R in (2.0e4, 2.0e5)
            E = scattered_field(p, q, λ, R .* n̂s)
            push!(e, abs(dot3(θ̂s, E) / predθ(R) - 1) + abs(dot3(φ̂s, E) / predφ(R) - 1))
        end
        @test e[1] < 2e-3              # close already at R = 2e4
        @test e[2] < e[1] / 5          # and converging as O(1/R)
    end
end
