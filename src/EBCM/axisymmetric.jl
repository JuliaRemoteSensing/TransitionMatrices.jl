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
transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng; stable = false) where {T, CT}
```

Calculate the T-Matrix for a given scatterer and wavelength.

Parameters:

- `s`: the axisymmetric scatterer.
- `λ`: the wavelength.
- `nₘₐₓ`: the maximum order of the T-Matrix.
- `Ng`: the number of Gauss-Legendre quadrature points to be used.
- `stable`: when `true`, assemble the `𝐔`-matrix integrals with the
  cancellation-free `F⁺` formulation of Somerville, Auguié & Le Ru, JQSRT 123
  (2013). The standard integrands lose all precision for spheroids of high aspect
  ratio (the irregular `χ_n·ψ_k` products develop huge Laurent terms that should
  cancel on integration but do not numerically); `stable=true` removes that
  cancellation, recovering a relative accuracy of about `1e-9` in `Float64`
  regardless of aspect ratio. It is **only valid for `Spheroid`** (the
  cancellation relies on the spheroid surface) and costs roughly 2–3× the default
  assembly, so it is opt-in. For it to help, `nₘₐₓ` must be large enough that
  `nₘₐₓ+1 ≳ k·c + 15`, where `c` is the largest semi-axis — comparable to the
  order needed for convergence at that size anyway. The remaining `Float64`
  round-off (in particular a `~1e-9` floor as the refractive index `s → 1`) is
  orthogonal to the cancellation and is lowered by using an extended-precision
  element type: `stable=true` with `Double64` reaches `~1e-25` at high aspect
  ratio.

Returns:

- `𝐓`: an `AxisymmetricTransitionMatrix` struct representing the T-Matrix.
"""
function transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng;
        zerofn = () -> zero(CT), stable = false) where {T, CT}
    stable && !(s isa Spheroid) &&
        throw(ArgumentError("stable=true is only valid for spheroids: the \
              cancellation-free F⁺ integrands rely on the spheroid surface making \
              the divergent Laurent terms integrate to zero (Somerville et al. 2013)."))
    𝐓 = Vector{Matrix{CT}}(undef, nₘₐₓ + 1)
    𝐓[1],
    cache = transition_matrix_m₀(s, λ, nₘₐₓ, Ng; zerofn = zerofn, reuse = true,
        stable = stable)
    for m in 1:nₘₐₓ
        𝐓[m + 1] = transition_matrix_m(m, s, λ, nₘₐₓ, Ng; zerofn = zerofn, cache = cache,
            stable = stable)
    end

    AxisymmetricTransitionMatrix{CT, nₘₐₓ, typeof(𝐓), T}(𝐓)
end

# Factor 𝐐 once so the factorization can be reused for the value block and
# every Jacobian slice (𝐐 depends only on 𝐏 and 𝐔, not on the differentiated
# parameter). Arblib matrices have no generic `lu`/`\`, so fall back to their
# dedicated `inv` and right-multiply. The `_ebcm_rdiv` methods dispatch on
# whether the second argument is an explicit inverse (a `Matrix`) or an `lu`
# factorization (which is not an `AbstractMatrix`).
_ebcm_factor(𝐐::AbstractMatrix{<:Union{Arb, Acb}}) = inv(𝐐)
_ebcm_factor(𝐐::AbstractMatrix) = lu(𝐐)

_ebcm_rdiv(𝐀, 𝐐⁻¹::AbstractMatrix) = 𝐀 * 𝐐⁻¹
_ebcm_rdiv(𝐀, F) = 𝐀 / F

function 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)
    𝐐 = @. 𝐏 + 1im * 𝐔
    # 𝐓 = -𝐏 𝐐⁻¹, via a factorization instead of an explicit inverse.
    return -_ebcm_rdiv(𝐏, _ebcm_factor(𝐐))
end

function ∂𝐓_from_𝐏_and_𝐔(𝐏, 𝐔, ∂𝐏, ∂𝐔)
    𝐐 = @. 𝐏 + 1im * 𝐔
    ∂𝐐 = @. ∂𝐏 + 1im * ∂𝐔
    F = _ebcm_factor(𝐐)
    𝐓 = -_ebcm_rdiv(𝐏, F)
    # ∂𝐓 = -∂𝐏 𝐐⁻¹ + 𝐏 𝐐⁻¹ ∂𝐐 𝐐⁻¹ = -(∂𝐏 + 𝐓 ∂𝐐) 𝐐⁻¹  (since 𝐏 𝐐⁻¹ = -𝐓),
    # reusing the single factorization `F` and the already-computed `𝐓`.
    return -_ebcm_rdiv(∂𝐏 + 𝐓 * ∂𝐐, F)
end

function _axisymmetric_transition_matrix_from_blocks(𝐓s::AbstractVector)
    N = length(𝐓s) - 1
    CT = eltype(first(𝐓s))
    T = real(CT)
    AxisymmetricTransitionMatrix{CT, N, typeof(𝐓s), T}(𝐓s)
end

function ebcm_transition_matrix_from_matrices(𝐏s::AbstractVector, 𝐔s::AbstractVector)
    length(𝐏s) == length(𝐔s) ||
        throw(ArgumentError("EBCM P and U block counts must match"))
    𝐓s = [𝐓_from_𝐏_and_𝐔(𝐏, 𝐔) for (𝐏, 𝐔) in zip(𝐏s, 𝐔s)]
    _axisymmetric_transition_matrix_from_blocks(𝐓s)
end

function ∂ebcm_transition_matrix_from_matrices(𝐏s::AbstractVector,
        𝐔s::AbstractVector,
        ∂𝐏s::AbstractVector,
        ∂𝐔s::AbstractVector)
    length(𝐏s) == length(𝐔s) == length(∂𝐏s) == length(∂𝐔s) ||
        throw(ArgumentError("EBCM P, U, ∂P, and ∂U block counts must match"))
    ∂𝐓s = [∂𝐓_from_𝐏_and_𝐔(𝐏, 𝐔, ∂𝐏, ∂𝐔)
           for (𝐏, 𝐔, ∂𝐏, ∂𝐔) in zip(𝐏s, 𝐔s, ∂𝐏s, ∂𝐔s)]
    _axisymmetric_transition_matrix_from_blocks(∂𝐓s)
end

"""
```
ebcm_matrices_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `P` and `U` matrices for the `m=0` EBCM block.
"""
function ebcm_matrices_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
        Ng; zerofn = () -> zero(CT), reuse = false,
        stable = false) where {T, CT}
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

    # For spheroids the off-diagonal `𝐔` integrands lose all precision at high
    # aspect ratio (the χ_n·ψ_k products have huge cancelling negative powers).
    # When `stable`, replace those products by the cancellation-free `F⁺` matrix
    # (Somerville et al. 2013); the diagonal (n=n′) has no cancellation and `𝐏`
    # is regular, so both keep the standard direct products. The x-independent
    # last-row series coefficients are computed once and shared across all points.
    Fmats = if stable
        lastrow = _F⁺_lastrow(s.m, nₘₐₓ, k * rₘₐₓ)
        [_F⁺_matrix(s.m, k * r[i], nₘₐₓ; lastrow) for i in 1:ng]
    else
        nothing
    end

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

                if stable
                    F = Fmats[i]
                    xi = k * r[i]
                    χψ = _Fp(F, n, n′) / xi
                    UL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * χψ
                    UL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * χψ
                    UL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * (_L⁷⁺_mat(F, n, n′) / xi)
                    UL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * (_L⁸⁺_mat(F, n, n′) / xi)
                else
                    UL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * χ[i, n] * ψₛ[i, n′]
                    UL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * χ[i, n] * ψₛ[i, n′]
                    UL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                           (χ′[i, n] * ψₛ′[i, n′] +
                            a[n] * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                    UL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                           (χ′[i, n] * ψₛ′[i, n′] +
                            a[n′] * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                end
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

    if reuse
        # `Fmats` (the per-quadrature-point F⁺ matrices) depend only on s.m, k and
        # r[i] — not on m — so they are computed once here and reused for all m>0.
        cache = x, w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, χₛ, χₛ′, Fmats
        return 𝐏, 𝐔, cache
    end

    return 𝐏, 𝐔
end

"""
```
transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `m=0` block of the T-Matrix for a given axisymmetric scatterer.
"""
function transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
        Ng; zerofn = () -> zero(CT), reuse = false,
        stable = false) where {T, CT}
    if reuse
        𝐏, 𝐔, cache = ebcm_matrices_m₀(s, λ, nₘₐₓ, Ng; zerofn, reuse, stable)
        return 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔), cache
    end

    𝐏, 𝐔 = ebcm_matrices_m₀(s, λ, nₘₐₓ, Ng; zerofn, reuse, stable)
    return 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)
end

"""
```
ebcm_matrices_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `P` and `U` matrices for the `m`-th EBCM block.
"""
function ebcm_matrices_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
        Ng; zerofn = () -> zero(CT), cache = nothing,
        stable = false) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    k = 2 * T(π) / λ
    nₘᵢₙ = max(1, m)
    nn = nₘₐₓ - nₘᵢₙ + 1
    sym = has_symmetric_plane(s)
    ng = sym ? Ng ÷ 2 : Ng

    if !isnothing(cache)
        _, w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, _, _, Fmats = cache
        return _transition_matrix_m_core(m, s, k, nₘₐₓ, Ng, nₘᵢₙ, nn, ng, sym, zerofn,
            w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, Fmats)
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

        Threads.@threads for i in 1:ng
            kₛr = k * s.m * r[i]
            ricattibesselj!(view(ψₛ, i, :), view(ψₛ′, i, :), view(zₛ, :, i), nₘₐₓ, nₑₓₜᵣₐ,
                kₛr)
        end

        Fmats = if stable
            lastrow = _F⁺_lastrow(s.m, nₘₐₓ, k * rₘₐₓ)
            [_F⁺_matrix(s.m, k * r[i], nₘₐₓ; lastrow) for i in 1:ng]
        else
            nothing
        end

        return _transition_matrix_m_core(m, s, k, nₘₐₓ, Ng, nₘᵢₙ, nn, ng, sym, zerofn,
            w, r, r′, ϑ, a, A, ψ, ψ′, χ, χ′, ψₛ, ψₛ′, Fmats)
    end
end

function _transition_matrix_m_core(m, s::AbstractAxisymmetricShape{T, CT}, k, nₘₐₓ, Ng,
        nₘᵢₙ, nn, ng, sym, zerofn, w, r, r′, ϑ, a, A, ψ,
        ψ′, χ, χ′, ψₛ, ψₛ′, Fmats = nothing) where {T, CT}
    stable = Fmats !== nothing
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

                if stable
                    F = Fmats[i]
                    xi = k * r[i]
                    # Cancellation-free K¹/K² integrands (Somerville, Auguié & Le Ru,
                    # JQSRT 123 (2013), Eqs. 53–54): K¹ ← [x·χ_n·ψ′_{n′}]⁺/x (their
                    # Eq. 59), K² ← [x·χ′_n·ψ_{n′}]⁺/x (their Eq. 60).
                    UK₁ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] *
                           (_xχψ′⁺_mat(F, n, n′) / xi)
                    UK₂ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] *
                           (_xχ′ψ⁺_mat(F, n, n′) / xi)
                else
                    UK₁ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] * χ[i, n] * ψₛ′[i, n′]
                    UK₂ += w[i] * k * r′[i] * 𝜋[i, n] * d[i, n′] * χ′[i, n] * ψₛ[i, n′]
                end
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

                    if stable
                        F = Fmats[i]
                        xi = k * r[i]
                        χψ = _Fp(F, n, n′) / xi
                        UL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * χψ
                        UL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * χψ
                        UL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                               (_L⁷⁺_mat(F, n, n′) / xi)
                        UL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                               (_L⁸⁺_mat(F, n, n′) / xi)
                    else
                        UL₁ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] * χ[i, n] * ψₛ[i, n′]
                        UL₂ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] * χ[i, n] * ψₛ[i, n′]
                        UL₇ += w[i] * k * r′[i] * τ[i, n] * d[i, n′] *
                               (χ′[i, n] * ψₛ′[i, n′] +
                                n * (n + 1) * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                        UL₈ += w[i] * k * r′[i] * d[i, n] * τ[i, n′] *
                               (χ′[i, n] * ψₛ′[i, n′] +
                                n′ * (n′ + 1) * χ[i, n] * ψₛ[i, n′] / (s.m * (k * r[i])^2))
                    end
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

    return 𝐏, 𝐔
end

"""
```
transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}
```

Calculate the `m`-th block of the T-Matrix for a given axisymmetric scatterer.
"""
function transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
        Ng; zerofn = () -> zero(CT), cache = nothing,
        stable = false) where {T, CT}
    𝐏, 𝐔 = ebcm_matrices_m(m, s, λ, nₘₐₓ, Ng; zerofn, cache, stable)
    return 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)
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
