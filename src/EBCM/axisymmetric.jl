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

Calculate the T-Matrix for a given scatterer and wavelength, using the given maximum order `nₘₐₓ` and number of Gauss-Legendre quadrature points `Ng`.

Parameters:

- `s`: the axisymmetricsScatterer.
- `λ`: the wavelength.
- `nₘₐₓ`: the maximum order of the T-Matrix.
- `Ng`: the number of Gauss-Legendre quadrature points to be used.

Returns:

- `𝐓`: an `AxisymmetricTransitionMatrix` struct representing the T-Matrix.
"""
function transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng;
                           zerofn = () -> zero(CT)) where {T, CT}
    𝐓 = Vector{Matrix{CT}}(undef, nₘₐₓ + 1)
    𝐓[1], cache = transition_matrix_m₀(s, λ, nₘₐₓ, Ng; zerofn = zerofn, reuse = true)
    for m in 1:nₘₐₓ
        𝐓[m + 1] = transition_matrix_m(m, s, λ, nₘₐₓ, Ng; zerofn = zerofn, cache = cache)
    end

    AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(𝐓), T}(𝐓)
end

function 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)
    𝐐 = @. 𝐏 + 1im * 𝐔
    𝐓 = -𝐏 * inv(𝐐)
end

## Hand-written version autodiff, slower than ForwardDiff.jl
# function 𝐓_from_𝐏_and_𝐔(𝐏::Matrix{Complex{ForwardDiff.Dual{ForwardDiff.Tag{F, T}, T, N}}},
#                         𝐔) where {F, T, N}
#     𝐐 = @. 𝐏 + 1im * 𝐔
#     𝐐r = real.(𝐐)
#     𝐐i = imag.(𝐐)
#     Qv = @. complex(ForwardDiff.value(𝐐r), ForwardDiff.value(𝐐i))
#     Qv⁻¹ = inv(Qv)
#     𝐏r = real.(𝐏)
#     𝐏i = imag.(𝐏)
#     Pv = @. complex(ForwardDiff.value(𝐏r), ForwardDiff.value(𝐏i))
#     Tv = -Pv * Qv⁻¹
#     ∂Pr = ForwardDiff.partials.(𝐏r)
#     ∂Pi = ForwardDiff.partials.(𝐏i)
#     ∂Qr = ForwardDiff.partials.(𝐐r)
#     ∂Qi = ForwardDiff.partials.(𝐐i)

#     ∂P = [map(zip(∂Pr, ∂Pi)) do (r, i)
#               complex(r[j], i[j])
#           end
#           for j in 1:N]
#     ∂Q = [map(zip(∂Qr, ∂Qi)) do (r, i)
#               complex(r[j], i[j])
#           end
#           for j in 1:N]

#     ∂T = map(zip(∂P, ∂Q)) do (∂Pj, ∂Qj)
#         -(∂Pj + Tv * ∂Qj) * Qv⁻¹
#     end

#     DT = eltype(𝐐r)
#     𝐓 = similar(𝐐)
#     map!(𝐓, CartesianIndices(Tv)) do ij
#         complex(DT(real(Tv[ij]),
#                    ForwardDiff.Partials(tuple([real(∂T[i][ij]) for i in 1:N]...))),
#                 DT(imag(Tv[ij]),
#                    ForwardDiff.Partials(tuple([imag(∂T[i][ij]) for i in 1:N]...))))
#     end
# end

"""
```
transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `m=0` block of the T-Matrix for a given axisymmetric scatterer.
"""
function transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
                              Ng; zerofn = () -> zero(CT), reuse = false) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    k = 2 * T(π) / λ
    sym = has_symmetric_plane(s)
    ng = sym ? Ng ÷ 2 : Ng

    x, w, r, r′ = gaussquad(s, Ng)
    ϑ = acos.(x)

    a = [n * (n + 1) for n in 1:nₘₐₓ]
    A = [√(T(2n + 1) / (2n * (n + 1))) for n in 1:nₘₐₓ]
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
        ricattibesselj!(view(ψₛ, i, :), view(ψₛ′, i, :), view(zₛ, :, i), nₘₐₓ, nₑₓₜᵣₐ,
                        kₛr)
        ricattibessely!(view(χₛ, i, :), view(χₛ′, i, :), nₘₐₓ, kₛr)
    end

    d = OffsetArray(zeros(T, Ng, nₘₐₓ + 1), 1:Ng, 0:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    Threads.@threads for i in eachindex(ϑ)
        wigner_d_recursion!(view(d, i, :), 0, 0, nₘₐₓ, ϑ[i];
                            deriv = view(τ, i, :))

        for n in 0:nₘₐₓ
            𝜋[i, n] = pi_func(T, 0, n, ϑ[i]; d = d[i, n])
        end
    end

    𝐏 = zeros(CT, 2nₘₐₓ, 2nₘₐₓ)
    𝐏₁₁ = view(𝐏, 1:nₘₐₓ, 1:nₘₐₓ)
    𝐏₂₂ = view(𝐏, (nₘₐₓ + 1):(2nₘₐₓ), (nₘₐₓ + 1):(2nₘₐₓ))
    𝐔 = zeros(CT, 2nₘₐₓ, 2nₘₐₓ)
    𝐔₁₁ = view(𝐔, 1:nₘₐₓ, 1:nₘₐₓ)
    𝐔₂₂ = view(𝐔, (nₘₐₓ + 1):(2nₘₐₓ), (nₘₐₓ + 1):(2nₘₐₓ))

    Threads.@threads for (n, n′) in collect(Iterators.product(1:nₘₐₓ, 1:nₘₐₓ))
        if sym && isodd(n + n′)
            continue
        end

        if n != n′
            PL₁ = zerofn()
            PL₂ = zerofn()
            PL₇ = zerofn()
            PL₈ = zerofn()

            UL₁ = zerofn()
            UL₂ = zerofn()
            UL₇ = zerofn()
            UL₈ = zerofn()

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
            PL̃₁ = zerofn()
            PL̃₂ = zerofn()
            PL̃₃ = zerofn()

            UL̃₁ = zerofn()
            UL̃₂ = zerofn()
            UL̃₃ = zerofn()

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

    𝐓 = 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)

    if reuse
        cache = x, w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, χₛ, χₛ′
        return 𝐓, cache
    end

    return 𝐓
end

"""
```
transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `m`-th block of the T-Matrix for a given axisymmetric scatterer.
"""
function transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
                             Ng; zerofn = () -> zero(CT), cache = nothing) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    k = 2 * T(π) / λ
    nₘᵢₙ = max(1, m)
    nn = nₘₐₓ - nₘᵢₙ + 1
    sym = has_symmetric_plane(s)
    ng = sym ? Ng ÷ 2 : Ng

    if !isnothing(cache)
        x, w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, χₛ, χₛ′ = cache
    else
        x, w, r, r′ = gaussquad(s, Ng)
        ϑ = acos.(x)
        a = OffsetArray([T(n * (n + 1)) for n in nₘᵢₙ:nₘₐₓ], nₘᵢₙ:nₘₐₓ)
        A = OffsetArray([√(T(2n + 1) / (2n * (n + 1))) for n in nₘᵢₙ:nₘₐₓ], nₘᵢₙ:nₘₐₓ)

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
            ricattibesselj!(view(ψₛ, i, :), view(ψₛ′, i, :), view(zₛ, :, i), nₘₐₓ, nₑₓₜᵣₐ,
                            kₛr)
            ricattibessely!(view(χₛ, i, :), view(χₛ′, i, :), nₘₐₓ, kₛr)
        end
    end

    d = OffsetArray(zeros(T, Ng, nₘₐₓ - m + 1), 1:Ng, m:nₘₐₓ)
    𝜋 = similar(d)
    τ = similar(d)

    Threads.@threads for i in eachindex(ϑ)
        wigner_d_recursion!(view(d, i, :), 0, m, nₘₐₓ, ϑ[i];
                            deriv = view(τ, i, :))

        for n in nₘᵢₙ:nₘₐₓ
            𝜋[i, n] = pi_func(T, m, n, ϑ[i]; d = d[i, n])
        end
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
            PK₁ = zerofn()
            PK₂ = zerofn()

            UK₁ = zerofn()
            UK₂ = zerofn()

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
                PL₁ = zerofn()
                PL₂ = zerofn()
                PL₇ = zerofn()
                PL₈ = zerofn()

                UL₁ = zerofn()
                UL₂ = zerofn()
                UL₇ = zerofn()
                UL₈ = zerofn()

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
                PL̃₁ = zerofn()
                PL̃₂ = zerofn()
                PL̃₃ = zerofn()

                UL̃₁ = zerofn()
                UL̃₂ = zerofn()
                UL̃₃ = zerofn()

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

    𝐓 = 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)

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

function extinction_efficiency_m₀(T₀)
    nₘₐₓ = size(T₀, 1) ÷ 2
    Qᵉˣᵗ = sum((2n + 1) * real(T₀[n, n] + T₀[n + nₘₐₓ, n + nₘₐₓ]) for n in 1:nₘₐₓ)
    return Qᵉˣᵗ
end

function scattering_efficiency_m₀(T₀)
    nₘₐₓ = size(T₀, 1) ÷ 2
    Qˢᶜᵃ = sum((2n + 1) *
               real(T₀[n, n] * T₀[n, n]' + T₀[n + nₘₐₓ, n + nₘₐₓ] * T₀[n + nₘₐₓ, n + nₘₐₓ]')
               for n in 1:nₘₐₓ)
    return Qˢᶜᵃ
end

@doc raw"""
```
expansion_coefficients(𝐓, λ)
```

Calculate the expansion coefficients from a given T-Matrix.

Parameters:

- `𝐓`: The precalculated T-Matrix of a scatterer.
- `λ`: The wavelength.
"""
function expansion_coefficients(𝐓::AxisymmetricTransitionMatrix{CT, N, V, T},
                                λ) where {CT, N, V, T}
    Cˢᶜᵃ = Float64(scattering_cross_section(𝐓, λ))
    λ = Float64(λ)

    ci = OffsetArray([(1.0im)^(i % 4) for i in (-N):N], (-N):N)
    s = OffsetArray([Float64(2i + 1) for i in 0:(2N)], 0:(2N))
    ss = sqrt.(s)
    sig = OffsetArray([1 - 2 * (i % 2) for i in 0:(4N)], 0:(4N))

    T1 = OffsetArray(zeros(ComplexF64, 2N + 1, N), (-N):N, 1:N)
    T2 = OffsetArray(zeros(ComplexF64, 2N + 1, N), (-N):N, 1:N)
    A1 = zeros(ComplexF64, N)
    A2 = zeros(ComplexF64, N)
    B1 = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N), 0:(2N), (-N):N, 1:N)
    B2 = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N), 0:(2N), (-N):N, 1:N)

    wig_table_init(4N, 3)
    wig_temp_init(4N)

    for n in 1:N
        # Calculate T1 and T2
        for n′ in 1:N
            for m in 0:min(n, n′)
                T11 = 𝐓[m, n, m, n′, 1, 1]
                T12 = 𝐓[m, n, m, n′, 1, 2]
                T21 = 𝐓[m, n, m, n′, 2, 1]
                T22 = 𝐓[m, n, m, n′, 2, 2]
                T1[m, n′] = T11 + T12 + T21 + T22
                T2[m, n′] = T11 + T12 - T21 - T22

                if m != 0
                    T1[-m, n′] = T11 - T12 - T21 + T22
                    T2[-m, n′] = T11 - T12 + T21 - T22
                end
            end
        end

        for n₁ in 0:(N + n)
            # Calculate A1 and A2
            for n′ in max(1, abs(n - n₁)):min(N, n₁ + n)
                A1[n′] = complex(0.0)
                A2[n′] = complex(0.0)
                for m₁ in (-min(n, n′)):min(n, n′)
                    cg = clebschgordan(n, m₁, n₁, 0, n′)
                    A1[n′] += cg * T1[m₁, n′]
                    A2[n′] += cg * T2[m₁, n′]
                end
                A1[n′] *= ci[n′ - n] / ss[n′]
                A2[n′] *= ci[n′ - n] / ss[n′]
            end

            # Calculate B1 and B2
            for m in max(1 - n₁, -n):min(n₁ + 1, n)
                for n′ in max(1, abs(n - n₁)):min(N, n₁ + n)
                    cg = clebschgordan(n, m, n₁, 1 - m, n′)
                    B1[n₁, m, n] += cg * A1[n′]
                    B2[n₁, m, n] += cg * A2[n′]
                end
            end
        end
    end

    # Calculate D
    D₀₀ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₀₋₀ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₂₂ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₂₋₂ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)

    for n in 1:N
        for n′ in 1:N
            for m in (-min(n, n′)):min(n, n′)
                for n₁ in abs(m - 1):(min(n, n′) + N)
                    D₀₀[m, n′, n] += s[n₁] * real(B1[n₁, m, n] * B1[n₁, m, n′]')
                    D₀₋₀[m, n′, n] += s[n₁] * real(B2[n₁, m, n] * B2[n₁, m, n′]')
                end
            end

            for m in max(-n, -n′ + 2):min(n, n′ + 2)
                for n₁ in abs(m - 1):(min(n, n′) + N)
                    D₂₂[m, n′, n] += s[n₁] * real(B1[n₁, m, n] * B1[n₁, 2 - m, n′]')
                    D₂₋₂[m, n′, n] += s[n₁] * real(B2[n₁, m, n] * B2[n₁, 2 - m, n′]')
                    D₀₂[m, n′, n] += s[n₁] * B2[n₁, m, n] * B1[n₁, 2 - m, n′]'
                end
            end
        end
    end

    h_const = λ^2 / (Cˢᶜᵃ * 4 * π)
    h = OffsetArray([s[l] * h_const * ss[n] / ss[n′]
                     for l in 0:(2N), n in 1:N, n′ in 1:N],
                    0:(2N),
                    1:N,
                    1:N)

    # Calculate g
    g₀₀ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₀₋₀ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₂₂ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₂₋₂ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))

    for l in 0:(2N)
        for n in 1:N
            for n′ in max(1, abs(n - l)):min(N, n + l)
                cg1 = clebschgordan(n, 1, l, 0, n′)
                sm₀₀ = 0.0
                sm₀₋₀ = 0.0
                for m in (-min(n, n′)):min(n, n′)
                    cg = clebschgordan(n, m, l, 0, n′)
                    sm₀₀ += cg * D₀₀[m, n′, n]
                    sm₀₋₀ += cg * D₀₋₀[m, n′, n]
                end
                g₀₀[l] += h[l, n, n′] * cg1 * sm₀₀
                g₀₋₀[l] += h[l, n, n′] * cg1 * sig[n + n′ + l] * sm₀₋₀

                if l >= 2
                    cg2 = clebschgordan(n, -1, l, 2, n′)
                    sm₂₂ = 0.0
                    sm₂₋₂ = 0.0
                    sm₀₂ = complex(0.0)
                    for m in max(-n, -n′ + 2):min(n, n′ + 2)
                        cg = clebschgordan(n, -m, l, 2, n′)
                        sm₂₂ += cg * D₂₂[m, n′, n]
                        sm₂₋₂ += cg * D₂₋₂[m, n′, n]
                        sm₀₂ += cg * D₀₂[m, n′, n]
                    end
                    g₂₂[l] += h[l, n, n′] * cg2 * sm₂₂
                    g₂₋₂[l] += h[l, n, n′] * cg2 * sig[n + n′ + l] * sm₂₋₂
                    g₀₂[l] += -h[l, n, n′] * cg1 * sm₀₂
                end
            end
        end
    end

    α₁ = g₀₀ + g₀₋₀
    α₂ = g₂₂ + g₂₋₂
    α₃ = g₂₂ - g₂₋₂
    α₄ = g₀₀ - g₀₋₀
    β₁ = 2real.(g₀₂)
    β₂ = 2imag.(g₀₂)

    wig_temp_free()
    wig_table_free()

    return α₁, α₂, α₃, α₄, β₁, β₂
end

@doc raw"""
```
scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θs)
```

Calculate the scatterering matrix elements from the given expansion coefficients.

Parameters:

- `α₁`, `α₂`, `α₃`, `α₄`, `β₁`, `β₂`: The precalculated expansion coefficients.
- `θs`: The scattering angles to be evaluated in degrees.
"""
function scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θs::AbstractVector)
    lmax = length(α₁) - 1
    θs = deg2rad.(θs)
    Nθ = length(θs)

    F = zeros(Nθ, 6)
    Threads.@threads for i in eachindex(θs)
        θ = θs[i]
        d₀₀ = wigner_d_recursion(0, 0, lmax, θ)
        d₂₂ = wigner_d_recursion(2, 2, lmax, θ)
        d₂₋₂ = wigner_d_recursion(2, -2, lmax, θ)
        d₀₂ = wigner_d_recursion(0, 2, lmax, θ)

        F₁₁ = sum(α₁[l] * d₀₀[l] for l in 0:lmax)
        F₂₂₊₃₃ = sum((α₂[l] + α₃[l]) * d₂₂[l] for l in 2:lmax)
        F₂₂₋₃₃ = sum((α₂[l] - α₃[l]) * d₂₋₂[l] for l in 2:lmax)
        F₂₂ = (F₂₂₊₃₃ + F₂₂₋₃₃) / 2
        F₃₃ = F₂₂₊₃₃ - F₂₂
        F₄₄ = sum(α₄[l] * d₀₀[l] for l in 0:lmax)
        F₁₂ = -sum(β₁[l] * d₀₂[l] for l in 2:lmax)
        F₃₄ = -sum(β₂[l] * d₀₂[l] for l in 2:lmax)

        F[i, :] .= F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄
    end

    return F
end

@doc raw"""
```
scattering_matrix(𝐓, λ, θs)
```

Calculate expansion coefficients first and then calculate scatterering matrix elements.

Parameters:

- `𝐓`: The transition matrix.
- `λ`: The wavelength.
- `θs`: The scattering angles to be evaluated in degrees.
"""
function scattering_matrix(𝐓::AxisymmetricTransitionMatrix, λ, θs::AbstractVector)
    α₁, α₂, α₃, α₄, β₁, β₂ = expansion_coefficients(𝐓, λ)
    return scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θs)
end

@testitem "Can calculate scattering matrix" begin
    s = Spheroid(1.0, 2.0, complex(1.311))
    λ = 2π
    𝐓 = transition_matrix(s, λ)
    θs = collect(0:180)
    F = scattering_matrix(𝐓, λ, θs)

    @test size(F) == (181, 6)
end
