# TODO: Add signature
@doc raw"""

Calculate the T-Matrix for a given scatterer and wavelength.
"""
function transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ) where {T, CT}
end

function transition_matrix_m₀(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
                              Ng) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    x, w = gausslegendre(T, Ng)
    ϑ = acos.(x)
    r = similar(x)
    r′ = similar(x)
    k = 2 * T(π) / λ
    radius_and_deriv!(r, r′, s, x)

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

function transition_matrix_m(m, s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ,
                             Ng) where {T, CT}
    @assert iseven(Ng) "Ng must be even!"

    x, w = gausslegendre(T, Ng)
    ϑ = acos.(x)
    r = similar(x)
    r′ = similar(x)
    k = 2 * T(π) / λ
    radius_and_deriv!(r, r′, s, x)

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

    @testset "Chebyshev" begin
        params = Iterators.product((0.5, 1.0, 5.0), (-0.5, 0.1, 0.9), (2, 3, 8),
                                   (1.311, 1.5 + 0.01im))
        nₘₐₓ = 10
        Ng = 200
        λ = 2π
        @testset "r₀ = $r₀, ε = $ε, n = $n, m = $m" for (r₀, ε, n, m) in params
            s = Chebyshev{Float64, ComplexF64}(r₀, ε, n, m)
            𝐓 = transition_matrix_m(0, s, λ, nₘₐₓ, Ng)
            𝐓₀ = transition_matrix_m₀(s, λ, nₘₐₓ, Ng)
            @test all(𝐓 .≈ 𝐓₀)
        end
    end
end
