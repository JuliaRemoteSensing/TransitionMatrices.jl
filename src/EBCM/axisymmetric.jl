struct AxisymmetricTransitionMatrix{CT, N, V <: AbstractVector{<:AbstractMatrix{CT}}, T} <:
       AbstractTransitionMatrix{CT, N}
    𝐓::V
end

Base.@propagate_inbounds function Base.getindex(axi::AxisymmetricTransitionMatrix{CT, N, V},
                                                m::Integer, n::Integer, m′::Integer,
                                                n′::Integer, p::Integer,
                                                p′::Integer) where {CT, N, V}
    if m != m′ || abs(m) > min(n, n′)
        zero(CT)
    else
        mₐ = abs(m)
        sig = m >= 0 ? 1 : (-1)^(p + p′)
        nn = N - max(1, mₐ) + 1
        n₁ = (p - 1) * nn + n - max(1, mₐ) + 1
        n₂ = (p′ - 1) * nn + n′ - max(1, mₐ) + 1
        axi.𝐓[mₐ + 1][n₁, n₂] * sig
    end
end

@doc raw"""
```
amplitude_matrix(axi::AxisymmetricTransitionMatrix{CT, N, V, T}, ϑᵢ, φᵢ, ϑₛ, φₛ;
                          λ = 2π,
                          rot::Union{Nothing, Rotation{3}} = nothing) where {CT, N, V, T}
```

Calculate the amplitude matrix of an axisymmetric scatterer.

Parameters:

- `axi`: the T-Matrix of the scatterer.
- `ϑᵢ`: the zenith angle of the incident wave.
- `φᵢ`: the azimuth angle of the incident wave.
- `ϑₛ`: the zenith angle of the scattered wave.
- `φₛ`: the azimuth angle of the scattered wave.
- `λ`: the wavelength of the incident wave in the host medium. Default to 2π.
- `rot`: the rotation of the scatterer.

"""
function amplitude_matrix(axi::AxisymmetricTransitionMatrix{CT, N, V, T}, ϑᵢ, φᵢ, ϑₛ, φₛ;
                          λ = 2π,
                          rot::Union{Nothing, Rotation{3}} = nothing) where {CT, N, V, T}
    k₁ = 2π / λ
    if isnothing(rot)
        α, β = zero(T), zero(T)
    else
        zyz = RotZYZ(rot)

        # gamma is useless for axisymmetric shapes
        α, β = T(zyz.theta1), T(zyz.theta2)
    end

    ϑᵢ = T(ϑᵢ)
    φᵢ = T(φᵢ)
    ϑₛ = T(ϑₛ)
    φₛ = T(φₛ)

    cosα = cos(α)
    sinα = sin(α)
    cosβ = cos(β)
    sinβ = sin(β)
    cosϑᵢ = cos(ϑᵢ)
    sinϑᵢ = sin(ϑᵢ)
    cosφ = cos(φᵢ - α)
    sinφ = sin(φᵢ - α)
    cosϑ₁ = cosϑᵢ * cosβ + sinϑᵢ * sinβ * cosφ
    ϑ₁ = acos(cosϑ₁)
    cosφ₁ = sinϑᵢ * cosβ * cosφ - cosϑᵢ * sinβ
    sinφ₁ = sinϑᵢ * sinφ
    φ₁ = atan(sinφ₁, cosφ₁)

    cosϑₛ = cos(ϑₛ)
    sinϑₛ = sin(ϑₛ)
    cosφ = cos(φₛ - α)
    sinφ = sin(φₛ - α)
    cosϑ₂ = cosϑₛ * cosβ + sinϑₛ * sinβ * cosφ
    ϑ₂ = acos(cosϑ₂)
    cosφ₂ = sinϑₛ * cosβ * cosφ - cosϑₛ * sinβ
    sinφ₂ = sinϑₛ * sinφ
    φ₂ = atan(sinφ₂, cosφ₂)

    𝐁 = @SMatrix [cosα*cosβ sinα*cosβ -sinβ
                  -sinα cosα 0
                  cosα*sinβ sinα*sinβ cosβ]

    cosφᵢ = cos(φᵢ)
    sinφᵢ = sin(φᵢ)
    cosφₛ = cos(φₛ)
    sinφₛ = sin(φₛ)

    𝐋ᵢ = @SMatrix [cosϑᵢ*cosφᵢ -sinφᵢ
                   cosϑᵢ*sinφᵢ cosφᵢ
                   -sinϑᵢ 0]
    𝐋ₛ = @SMatrix [cosϑₛ*cosφₛ -sinφₛ
                   cosϑₛ*sinφₛ cosφₛ
                   -sinϑₛ 0]

    sinϑ₁ = sin(ϑ₁)
    cosφ₁ = cos(φ₁)
    sinφ₁ = sin(φ₁)
    sinϑ₂ = sin(ϑ₂)
    cosφ₂ = cos(φ₂)
    sinφ₂ = sin(φ₂)

    𝐏₁ = @SMatrix [cosϑ₁*cosφ₁ cosϑ₁*sinφ₁ -sinϑ₁
                   -sinφ₁ cosφ₁ 0]
    𝐏₂ = @SMatrix [cosϑ₂*cosφ₂ cosϑ₂*sinφ₂ -sinϑ₂
                   -sinφ₂ cosφ₂ 0]

    𝐑₁ = 𝐏₁ * (𝐁 * 𝐋ᵢ)
    𝐑₂ = inv(𝐏₂ * (𝐁 * 𝐋ₛ))

    S₁₁, S₁₂, S₂₁, S₂₂ = zero(CT), zero(CT), zero(CT), zero(CT)
    φ = φ₂ - φ₁

    coeff = ([1.0im^((n′ - n - 1) & 3) *
              √T((2n + 1) * (2n′ + 1) / (n * (n + 1) * n′ * (n′ + 1)))
              for n in 1:N, n′ in 1:N])

    for m in 0:N
        π₁, τ₁ = wigner_d_recursion(T, 0, m, N, ϑ₁, deriv = true)
        π₂, τ₂ = wigner_d_recursion(T, 0, m, N, ϑ₂, deriv = true)

        for n in m:N
            π₁[n] = pi_func(T, m, n, ϑ₁; d = π₁[n])
            π₂[n] = pi_func(T, m, n, ϑ₂; d = π₂[n])
        end

        nₘ = N - max(1, m) + 1
        offset = N - nₘ
        cosmφ = cos(m * φ)
        sinmφ = sin(m * φ)

        for n′ in 1:nₘ, n in 1:nₘ
            T₁₁ = axi.𝐓[m + 1][n, n′]
            T₂₂ = axi.𝐓[m + 1][n + nₘ, n′ + nₘ]
            if m == 0
                S₁₁ += coeff[n + offset, n′ + offset] * τ₂[n + offset] * τ₁[n′ + offset] *
                       T₂₂
                S₂₂ += coeff[n + offset, n′ + offset] * τ₂[n + offset] * τ₁[n′ + offset] *
                       T₁₁
            else
                T₁₂ = axi.𝐓[m + 1][n, n′ + nₘ]
                T₂₁ = axi.𝐓[m + 1][n + nₘ, n′]
                c₁ = coeff[n + offset, n′ + offset] * 2cosmφ
                c₂ = coeff[n + offset, n′ + offset] * 2sinmφ

                D₁₁ = π₂[n + offset] * π₁[n′ + offset]
                D₁₂ = π₂[n + offset] * τ₁[n′ + offset]
                D₂₁ = τ₂[n + offset] * π₁[n′ + offset]
                D₂₂ = τ₂[n + offset] * τ₁[n′ + offset]
                S₁₁ += c₁ * (T₁₁ * D₁₁ + T₁₂ * D₁₂ + T₂₁ * D₂₁ + T₂₂ * D₂₂)
                S₁₂ += c₂ * (T₁₁ * D₁₂ + T₁₂ * D₁₁ + T₂₁ * D₂₂ + T₂₂ * D₂₁)
                S₂₁ -= c₂ * (T₁₁ * D₂₁ + T₁₂ * D₂₂ + T₂₁ * D₁₁ + T₂₂ * D₁₂)
                S₂₂ += c₁ * (T₁₁ * D₂₂ + T₁₂ * D₂₁ + T₂₁ * D₁₂ + T₂₂ * D₁₁)
            end
        end
    end

    𝐒 = @SMatrix [S₁₁ S₁₂
                  S₂₁ S₂₂]

    return 𝐑₂ * (𝐒 * 𝐑₁) / k₁
end

@testitem "𝐒(axi, euler) ≡ 𝐒(rotate(axi, euler))" begin
    using Rotations: RotZYZ
    using TransitionMatrices: Spheroid, amplitude_matrix, rotate, transition_matrix

    @testset "Spheroid" begin
        params = Iterators.product((0.5, 1.0, 5.0), (1.0,), (1.311, 1.5 + 0.01im),
                                   (0.0, 0.5), (0.0, 0.5))

        @testset "a = $a, c = $c, m = $m, α = $α, β = $β" for (a, c, m, α, β) in params
            s = Spheroid{Float64, ComplexF64}(a, c, m)
            𝐓 = transition_matrix(s, 2π, 15, 200)
            𝐓r = rotate(𝐓, RotZYZ(α, β, 0))

            𝐒 = amplitude_matrix(𝐓, 0.0, 0.3, π / 2, 0.5; rot = RotZYZ(α, β, 0))
            𝐒r = amplitude_matrix(𝐓r, 0.0, 0.3, π / 2, 0.5)

            @test all(𝐒 .≈ 𝐒r)
        end
    end
end

@doc raw"""
```
transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the T-Matrix for a given scatterer and wavelength.

Parameters:

- `s`: the axisymmetricsScatterer.
- `λ`: the wavelength.
- `nₘₐₓ`: the maximum order of the T-Matrix.
- `Ng`: the number of Gauss-Legendre quadrature points to be used.

Returns:

- `𝐓`: an `AxisymmetricTransitionMatrix` struct representing the T-Matrix.
"""
function transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
    𝐓 = Vector{Matrix{CT}}(undef, nₘₐₓ + 1)
    𝐓[1] = transition_matrix_m₀(s, λ, nₘₐₓ, Ng)
    Threads.@threads for m in 1:nₘₐₓ
        𝐓[m + 1] = transition_matrix_m(m, s, λ, nₘₐₓ, Ng)
    end

    AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(𝐓), T}(𝐓)
end

"""
```
transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `m=0` block of the T-Matrix for a given axisymmetric scatterer.
"""
function transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
                              Ng) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    x, w, r, r′ = gaussquad(s, Ng)
    ϑ = acos.(x)
    k = 2 * T(π) / λ

    a = [n * (n + 1) for n in 1:nₘₐₓ]
    A = [√(T(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]
    d = OffsetArray(zeros(T, Ng, nₘₐₓ + 1), 1:Ng, 0:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    for i in eachindex(ϑ)
        wigner_d_recursion!(view(d, i, :), 0, 0, nₘₐₓ, ϑ[i];
                            deriv = view(τ, i, :))

        for n in 0:nₘₐₓ
            𝜋[i, n] = pi_func(T, 0, n, ϑ[i]; d = d[i, n])
        end
    end

    sym = has_symmetric_plane(s)
    ng = sym ? Ng ÷ 2 : Ng

    rₘₐₓ = maximum(r)
    nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, k * rₘₐₓ)
    ψ = zeros(T, ng, nₘₐₓ)
    z = zeros(T, nₘₐₓ + nₑₓₜᵣₐ, ng)
    ψ′ = similar(ψ)
    χ = similar(ψ)
    χ′ = similar(ψ)

    Threads.@threads for i in 1:ng
        kr = k * r[i]
        ricattibesselj!(view(ψ, i, :), view(ψ′, i, :), view(z, :, i), nₘₐₓ, nₑₓₜᵣₐ, kr)
        ricattibessely!(view(χ, i, :), view(χ′, i, :), nₘₐₓ, kr)
    end

    nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, s.m * k * rₘₐₓ)
    ψₛ = zeros(CT, ng, nₘₐₓ)
    zₛ = zeros(CT, nₘₐₓ + nₑₓₜᵣₐ, ng)
    ψₛ′ = similar(ψₛ)
    χₛ = similar(ψₛ)
    χₛ′ = similar(ψₛ)

    𝐏 = zeros(CT, 2nₘₐₓ, 2nₘₐₓ)
    𝐏₁₁ = view(𝐏, 1:nₘₐₓ, 1:nₘₐₓ)
    𝐏₂₂ = view(𝐏, (nₘₐₓ + 1):(2nₘₐₓ), (nₘₐₓ + 1):(2nₘₐₓ))
    𝐔 = zeros(CT, 2nₘₐₓ, 2nₘₐₓ)
    𝐔₁₁ = view(𝐔, 1:nₘₐₓ, 1:nₘₐₓ)
    𝐔₂₂ = view(𝐔, (nₘₐₓ + 1):(2nₘₐₓ), (nₘₐₓ + 1):(2nₘₐₓ))

    Threads.@threads for i in 1:ng
        kₛr = k * s.m * r[i]
        ricattibesselj!(view(ψₛ, i, :), view(ψₛ′, i, :), view(zₛ, :, i), nₘₐₓ, nₑₓₜᵣₐ, kₛr)
        ricattibessely!(view(χₛ, i, :), view(χₛ′, i, :), nₘₐₓ, kₛr)
    end

    Threads.@threads for (n, n′) in collect(Iterators.product(1:nₘₐₓ, 1:nₘₐₓ))
        if sym && isodd(n + n′)
            continue
        end

        if n != n′
            PL₁ = zero(CT)
            PL₂ = zero(CT)
            PL₇ = zero(CT)
            PL₈ = zero(CT)

            UL₁ = zero(CT)
            UL₂ = zero(CT)
            UL₇ = zero(CT)
            UL₈ = zero(CT)

            for i in 1:ng
                PL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * ψ[i, n] * ψₛ[i, n′]
                PL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * ψ[i, n] * ψₛ[i, n′]
                PL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                       (ψ′[i, n] * ψₛ′[i, n′] +
                        a[n] * ψ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                PL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                       (ψ′[i, n] * ψₛ′[i, n′] +
                        a[n′] * ψ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))

                UL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * χ[i, n] * ψₛ[i, n′]
                UL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * χ[i, n] * ψₛ[i, n′]
                UL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                       (χ′[i, n] * ψₛ′[i, n′] +
                        a[n] * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                UL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                       (χ′[i, n] * ψₛ′[i, n′] +
                        a[n′] * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
            end

            𝐏₁₁[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                         (a[n] * PL₂ - a[n′] * PL₁)
            𝐏₂₂[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                         (a[n] * PL₈ - a[n′] * PL₇)
            𝐔₁₁[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                         (a[n] * UL₂ - a[n′] * UL₁)
            𝐔₂₂[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                         (a[n] * UL₈ - a[n′] * UL₇)
        else
            PL̃₁ = zero(CT)
            PL̃₂ = zero(CT)
            PL̃₃ = zero(CT)

            UL̃₁ = zero(CT)
            UL̃₂ = zero(CT)
            UL̃₃ = zero(CT)

            for i in 1:ng
                PL̃₁ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                        (ψ′[i, n] * ψₛ[i, n] - s.m * ψ[i, n] * ψₛ′[i, n])
                PL̃₂ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                        (s.m * ψ′[i, n] * ψₛ[i, n] - ψ[i, n] * ψₛ′[i, n])
                PL̃₃ += w[i] * k * r′[i] * τ[i, n] * d[i, n] * ψ[i, n] * ψₛ[i, n] /
                        (s.m * (k * r[i])^2)

                UL̃₁ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                        (χ′[i, n] * ψₛ[i, n] - s.m * χ[i, n] * ψₛ′[i, n])
                UL̃₂ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                        (s.m * χ′[i, n] * ψₛ[i, n] - χ[i, n] * ψₛ′[i, n])
                UL̃₃ += w[i] * k * r′[i] * τ[i, n] * d[i, n] * χ[i, n] * ψₛ[i, n] /
                        (s.m * (k * r[i])^2)
            end

            𝐏₁₁[n, n] = -1im / s.m * A[n]^2 * PL̃₁
            𝐏₂₂[n, n] = -1im / s.m * A[n]^2 * (PL̃₂ + (s.m^2 - 1) * a[n] * PL̃₃)

            𝐔₁₁[n, n] = -1im / s.m * A[n]^2 * UL̃₁
            𝐔₂₂[n, n] = -1im / s.m * A[n]^2 * (UL̃₂ + (s.m^2 - 1) * a[n] * UL̃₃)
        end
    end

    𝐐 = @. 𝐏 + 1im * 𝐔
    𝐓 = -𝐏 * inv(𝐐)

    return 𝐓
end

"""
```
transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `m`-th block of the T-Matrix for a given axisymmetric scatterer.
"""
function transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
                             Ng) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    x, w, r, r′ = gaussquad(s, Ng)
    ϑ = acos.(x)
    k = 2 * T(π) / λ

    nₘᵢₙ = max(1, m)
    nn = nₘₐₓ - nₘᵢₙ + 1
    a = OffsetArray([n * (n + 1) for n in nₘᵢₙ:nₘₐₓ], nₘᵢₙ:nₘₐₓ)
    A = OffsetArray([√(T(2n + 1) / (2n * (n + 1))) for n in nₘᵢₙ:nₘₐₓ], nₘᵢₙ:nₘₐₓ)
    d = OffsetArray(zeros(T, Ng, nₘₐₓ - m + 1), 1:Ng, m:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    for i in eachindex(ϑ)
        wigner_d_recursion!(view(d, i, :), 0, m, nₘₐₓ, ϑ[i];
                            deriv = view(τ, i, :))

        for n in nₘᵢₙ:nₘₐₓ
            𝜋[i, n] = pi_func(T, m, n, ϑ[i]; d = d[i, n])
        end
    end

    sym = has_symmetric_plane(s)
    ng = sym ? Ng ÷ 2 : Ng

    rₘₐₓ = maximum(r)
    nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, k * rₘₐₓ)
    ψ = zeros(T, ng, nₘₐₓ)
    z = zeros(T, nₘₐₓ + nₑₓₜᵣₐ, ng)
    ψ′ = similar(ψ)
    χ = similar(ψ)
    χ′ = similar(ψ)

    Threads.@threads for i in 1:ng
        kr = k * r[i]
        ricattibesselj!(view(ψ, i, :), view(ψ′, i, :), view(z, :, i), nₘₐₓ, nₑₓₜᵣₐ, kr)
        ricattibessely!(view(χ, i, :), view(χ′, i, :), nₘₐₓ, kr)
    end

    nₑₓₜᵣₐ = estimate_ricattibesselj_extra_terms(nₘₐₓ, s.m * k * rₘₐₓ)
    ψₛ = zeros(CT, ng, nₘₐₓ)
    zₛ = zeros(CT, nₘₐₓ + nₑₓₜᵣₐ, ng)
    ψₛ′ = similar(ψₛ)
    χₛ = similar(ψₛ)
    χₛ′ = similar(ψₛ)

    Threads.@threads for i in 1:ng
        kₛr = k * s.m * r[i]
        ricattibesselj!(view(ψₛ, i, :), view(ψₛ′, i, :), view(zₛ, :, i), nₘₐₓ, nₑₓₜᵣₐ, kₛr)
        ricattibessely!(view(χₛ, i, :), view(χₛ′, i, :), nₘₐₓ, kₛr)
    end

    𝐏 = zeros(CT, 2nn, 2nn)
    𝐏₁₁ = OffsetArray(view(𝐏, 1:nn, 1:nn), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    𝐏₁₂ = OffsetArray(view(𝐏, 1:nn, (nn + 1):(2nn)), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    𝐏₂₁ = OffsetArray(view(𝐏, (nn + 1):(2nn), 1:nn), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    𝐏₂₂ = OffsetArray(view(𝐏, (nn + 1):(2nn), (nn + 1):(2nn)), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)

    𝐔 = zeros(CT, 2nn, 2nn)
    𝐔₁₁ = OffsetArray(view(𝐔, 1:nn, 1:nn), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    𝐔₁₂ = OffsetArray(view(𝐔, 1:nn, (nn + 1):(2nn)), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    𝐔₂₁ = OffsetArray(view(𝐔, (nn + 1):(2nn), 1:nn), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)
    𝐔₂₂ = OffsetArray(view(𝐔, (nn + 1):(2nn), (nn + 1):(2nn)), nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ)

    Threads.@threads for (n, n′) in collect(Iterators.product(nₘᵢₙ:nₘₐₓ, nₘᵢₙ:nₘₐₓ))
        if !(sym && iseven(n + n′))
            PK₁ = zero(CT)
            PK₂ = zero(CT)

            UK₁ = zero(CT)
            UK₂ = zero(CT)

            for i in 1:ng
                PK₁ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] * ψ[i, n] * ψₛ′[i, n′]
                PK₂ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] * ψ′[i, n] * ψₛ[i, n′]

                UK₁ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] * χ[i, n] * ψₛ′[i, n′]
                UK₂ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] * χ′[i, n] * ψₛ[i, n′]
            end

            𝐏₁₂[n, n′] = A[n] * A[n′] * (s.m^2 - 1) / s.m * PK₁
            𝐏₂₁[n, n′] = A[n] * A[n′] * (1 - s.m^2) / s.m * PK₂

            𝐔₁₂[n, n′] = A[n] * A[n′] * (s.m^2 - 1) / s.m * UK₁
            𝐔₂₁[n, n′] = A[n] * A[n′] * (1 - s.m^2) / s.m * UK₂
        end

        if !(sym && isodd(n + n′))
            if n != n′
                PL₁ = zero(CT)
                PL₂ = zero(CT)
                PL₇ = zero(CT)
                PL₈ = zero(CT)

                UL₁ = zero(CT)
                UL₂ = zero(CT)
                UL₇ = zero(CT)
                UL₈ = zero(CT)

                for i in 1:ng
                    PL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * ψ[i, n] * ψₛ[i, n′]
                    PL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * ψ[i, n] * ψₛ[i, n′]
                    PL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                           (ψ′[i, n] * ψₛ′[i, n′] +
                            n * (n + 1) * ψ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                    PL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                           (ψ′[i, n] * ψₛ′[i, n′] +
                            n′ * (n′ + 1) * ψ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))

                    UL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * χ[i, n] * ψₛ[i, n′]
                    UL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * χ[i, n] * ψₛ[i, n′]
                    UL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                           (χ′[i, n] * ψₛ′[i, n′] +
                            n * (n + 1) * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                    UL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                           (χ′[i, n] * ψₛ′[i, n′] +
                            n′ * (n′ + 1) * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                end

                𝐏₁₁[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                             (a[n] * PL₂ - a[n′] * PL₁)
                𝐏₂₂[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                             (a[n] * PL₈ - a[n′] * PL₇)

                𝐔₁₁[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                             (a[n] * UL₂ - a[n′] * UL₁)
                𝐔₂₂[n, n′] = 1im * A[n] * A[n′] * (s.m^2 - 1) / (s.m * (a[n] - a[n′])) *
                             (a[n] * UL₈ - a[n′] * UL₇)
            else
                PL̃₁ = zero(CT)
                PL̃₂ = zero(CT)
                PL̃₃ = zero(CT)

                UL̃₁ = zero(CT)
                UL̃₂ = zero(CT)
                UL̃₃ = zero(CT)

                for i in 1:ng
                    PL̃₁ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                            (ψ′[i, n] * ψₛ[i, n] - s.m * ψ[i, n] * ψₛ′[i, n])
                    PL̃₂ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                            (s.m * ψ′[i, n] * ψₛ[i, n] - ψ[i, n] * ψₛ′[i, n])
                    PL̃₃ += w[i] * k * r′[i] * τ[i, n] * d[i, n] * ψ[i, n] * ψₛ[i, n] /
                            (s.m * (k * r[i])^2)

                    UL̃₁ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                            (χ′[i, n] * ψₛ[i, n] - s.m * χ[i, n] * ψₛ′[i, n])
                    UL̃₂ += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) *
                            (s.m * χ′[i, n] * ψₛ[i, n] - χ[i, n] * ψₛ′[i, n])
                    UL̃₃ += w[i] * k * r′[i] * τ[i, n] * d[i, n] * χ[i, n] * ψₛ[i, n] /
                            (s.m * (k * r[i])^2)
                end

                𝐏₁₁[n, n] = -1im / s.m * A[n]^2 * PL̃₁
                𝐏₂₂[n, n] = -1im / s.m * A[n]^2 * (PL̃₂ + (s.m^2 - 1) * a[n] * PL̃₃)

                𝐔₁₁[n, n] = -1im / s.m * A[n]^2 * UL̃₁
                𝐔₂₂[n, n] = -1im / s.m * A[n]^2 * (UL̃₂ + (s.m^2 - 1) * a[n] * UL̃₃)
            end
        end
    end

    𝐐 = @. 𝐏 + 1im * 𝐔
    𝐓 = -𝐏 * inv(𝐐)

    return 𝐓
end

@testitem "transition_matrix_m should be equivalent to transition_matrix_m₀ when m = 0" begin
    using TransitionMatrices: Spheroid, Chebyshev, transition_matrix_m, transition_matrix_m₀

    @testset "Spheroid" begin
        params = Iterators.product((1.0, 2.0, 5.0), (0.9, 1.8, 4.5), (1.311, 1.5 + 0.01im))
        nₘₐₓ = 10
        Ng = 200
        λ = 2π
        @testset "a = $a, c = $c, m = $m" for (a, c, m) in params
            s = Spheroid{Float64, ComplexF64}(a, c, m)
            𝐓 = transition_matrix_m(0, s, λ, nₘₐₓ, Ng)
            𝐓₀ = transition_matrix_m₀(s, λ, nₘₐₓ, Ng)
            @test all(𝐓 .≈ 𝐓₀)
        end
    end

    @testset "Cylinder" begin
        params = Iterators.product((1.0, 2.0), (0.5, 2.0), (1.311, 1.5 + 0.01im))
        nₘₐₓ = 10
        Ng = 200
        λ = 2π
        @testset "r = $r, h = $h, m = $m" for (r, h, m) in params
            c = Cylinder{Float64, ComplexF64}(r, h, m)
            𝐓 = transition_matrix_m(0, c, λ, nₘₐₓ, Ng)
            𝐓₀ = transition_matrix_m₀(c, λ, nₘₐₓ, Ng)
            @test all(𝐓 .≈ 𝐓₀)
        end
    end

    @testset "Chebyshev" begin
        params = Iterators.product((0.5, 1.0, 5.0), (-0.5, 0.1, 0.9), (2, 3, 8),
                                   (1.311, 1.5 + 0.01im))
        nₘₐₓ = 10
        Ng = 200
        λ = 2π
        @testset "r₀ = $r₀, ε = $ε, n = $n, m = $m" for (r₀, ε, n, m) in params
            c = Chebyshev{Float64, ComplexF64}(r₀, ε, n, m)
            𝐓 = transition_matrix_m(0, c, λ, nₘₐₓ, Ng)
            𝐓₀ = transition_matrix_m₀(c, λ, nₘₐₓ, Ng)
            @test all(𝐓 .≈ 𝐓₀)
        end
    end
end

@doc raw"""
```
scattering_cross_section(axi::AxisymmetricTransitionMatrix{CT, N}, λ=2π) where {CT, N}
```

Calculate the scattering cross section per particle averaged over the uniform orientation distribution, according to Eq. (5.141) in Mishchenko et al. (2002).

```math
\left\langle C_{\text {sca }}\right\rangle=\frac{2 \pi}{k_1^2} \sum_{n=1}^{\infty} \sum_{n^{\prime}=1}^{\infty} \sum_{m=0}^{\min \left(n, n^{\prime}\right)} \sum_{k=1}^2 \sum_{l=1}^2\left(2-\delta_{m 0}\right)\left|T_{m n m n^{\prime}}^{k l}(P)\right|^2
```

Parameters:

- `𝐓`: the T-Matrix of the scatterer.
- `λ`: the wavelength of the incident wave in the host medium. Default to 2π.
"""
function scattering_cross_section(𝐓::AxisymmetricTransitionMatrix{CT, N, V, T},
                                  λ = 2π) where {CT, N, V, T}
    Cˢᶜᵃ = zero(T)
    for m in 0:N
        for p′ in 1:2, p in 1:2
            for n′ in max(m, 1):N, n in max(m, 1):N
                if m == 0
                    Cˢᶜᵃ += abs2(𝐓[m, n, m, n′, p, p′])
                else
                    Cˢᶜᵃ += 2 * abs2(𝐓[m, n, m, n′, p, p′])
                end
            end
        end
    end

    Cˢᶜᵃ * λ^2 / 2π
end

@testitem "scattering cross section should be the same when calculating for axisymmetric scatterers using the general method" begin
    using TransitionMatrices: Spheroid, TransitionMatrix, transition_matrix,
                              scattering_cross_section

    s = Spheroid(1.0, 0.5, 1.5 + 0.01im)
    𝐓 = transition_matrix(s, 2π, 5, 40)
    Cˢᶜᵃ = scattering_cross_section(𝐓)
    Cˢᶜᵃ′ = scattering_cross_section(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓))
    @test Cˢᶜᵃ ≈ Cˢᶜᵃ′
end

@doc raw"""
```
extinction_cross_section(axi::AxisymmetricTransitionMatrix{CT, N}, λ=2π) where {CT, N}
```

Calculate the extinction cross section per particle averaged over the uniform orientation distribution, according to Eq. (5.107) in Mishchenko et al. (2002).

```math
\left\langle C_{\text {ext }}\right\rangle=-\frac{2 \pi}{k_1^2} \operatorname{Re} \sum_{n=1}^{\infty} \sum_{m=0}^n\left(2-\delta_{m 0}\right)\left[T_{m n m n}^{11}(P)+T_{m n m n}^{22}(P)\right]
```

Parameters:

- `𝐓`: the T-Matrix of the scatterer.
- `λ`: the wavelength of the incident wave in the host medium. Default to 2π.
"""
function extinction_cross_section(𝐓::AxisymmetricTransitionMatrix{CT, N, V, T},
                                  λ = 2π) where {CT, N, V, T}
    Cᵉˣᵗ = zero(CT)
    for m in 0:N
        coeff = m == 0 ? 1 : 2
        for n in max(m, 1):N
            Cᵉˣᵗ += coeff * (𝐓[m, n, m, n, 1, 1] + 𝐓[m, n, m, n, 2, 2])
        end
    end

    -real(Cᵉˣᵗ) * λ^2 / 2π
end

@testitem "extinction cross section should be the same when calculating for axisymmetric scatterers using the general method" begin
    using TransitionMatrices: Spheroid, TransitionMatrix, transition_matrix,
                              extinction_cross_section

    s = Spheroid(1.0, 0.5, 1.5 + 0.01im)
    𝐓 = transition_matrix(s, 2π, 5, 40)
    Cᵉˣᵗ = extinction_cross_section(𝐓)
    Cᵉˣᵗ′ = extinction_cross_section(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓))
    @test Cᵉˣᵗ ≈ Cᵉˣᵗ′
end
