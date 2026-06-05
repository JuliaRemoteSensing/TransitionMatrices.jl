struct AxisymmetricTransitionMatrix{CT, N, V <: AbstractVector{<:AbstractMatrix{CT}}, T} <:
       AbstractTransitionMatrix{CT, N}
    𝐓::V
end

Base.@propagate_inbounds function Base.getindex(
        axi::AxisymmetricTransitionMatrix{CT, N, V},
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
scattering_cross_section(axi::AxisymmetricTransitionMatrix{CT, N}, λ=2π) where {CT, N}
```

Calculate the scattering cross section per particle averaged over the uniform orientation distribution, according to Mishchenko et al. (2002), Eq. (5.141).

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

Calculate the extinction cross section per particle averaged over the uniform orientation distribution, according to Mishchenko et al. (2002), Eq. (5.107).

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
               real(T₀[n, n] * T₀[n, n]' +
                    T₀[n + nₘₐₓ, n + nₘₐₓ] * T₀[n + nₘₐₓ, n + nₘₐₓ]')
    for n in 1:nₘₐₓ)
    return Qˢᶜᵃ
end

@doc raw"""
```
expansion_coefficients(𝐓::AxisymmetricTransitionMatrix{CT, N, V, T}, λ) where {CT, N, V, T}
```

Calculate the expansion coefficients from an axisymmetric T-Matrix. Translated from Mishchenko et al.'s Fortran code.

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

    T1 = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    T2 = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)
    A1 = zeros(ComplexF64, N, N)
    A2 = zeros(ComplexF64, N, N)
    B1 = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N), 0:(2N), (-N):N, 1:N)
    B2 = OffsetArray(zeros(ComplexF64, 2N + 1, 2N + 1, N), 0:(2N), (-N):N, 1:N)

    tid_offset = VERSION >= v"1.11" ? Threads.nthreads(:interactive) : 0
    wig_table_init(4N, 3)
    @sync for i in 1:Threads.nthreads()
        StableTasks.@spawnat i + tid_offset wig_thread_temp_init(4N)
    end

    @debug "Calculating B..."
    Threads.@threads for n in 1:N
        @debug "n = $n..."

        # Calculate T1 and T2
        for n′ in 1:N
            for m in 0:min(n, n′)
                T11 = 𝐓[m, n, m, n′, 1, 1]
                T12 = 𝐓[m, n, m, n′, 1, 2]
                T21 = 𝐓[m, n, m, n′, 2, 1]
                T22 = 𝐓[m, n, m, n′, 2, 2]
                T1[m, n′, n] = T11 + T12 + T21 + T22
                T2[m, n′, n] = T11 + T12 - T21 - T22

                if m != 0
                    T1[-m, n′, n] = T11 - T12 - T21 + T22
                    T2[-m, n′, n] = T11 - T12 + T21 - T22
                end
            end
        end

        for n₁ in 0:(N + n)
            # Calculate A1 and A2
            for n′ in max(1, abs(n - n₁)):min(N, n₁ + n)
                a₁ = 0.0im
                a₂ = 0.0im
                for m₁ in (-min(n, n′)):min(n, n′)
                    cg = clebschgordan(n, m₁, n₁, 0, n′)
                    a₁ += cg * T1[m₁, n′, n]
                    a₂ += cg * T2[m₁, n′, n]
                end
                a₁ *= ci[n′ - n] / ss[n′]
                a₂ *= ci[n′ - n] / ss[n′]
                A1[n′, n] = a₁
                A2[n′, n] = a₂
            end

            # Calculate B1 and B2
            for m in max(1 - n₁, -n):min(n₁ + 1, n)
                b₁ = 0.0im
                b₂ = 0.0im
                for n′ in max(1, abs(n - n₁)):min(N, n₁ + n)
                    cg = clebschgordan(n, m, n₁, 1 - m, n′)
                    b₁ += cg * A1[n′, n]
                    b₂ += cg * A2[n′, n]
                end
                B1[n₁, m, n] = b₁
                B2[n₁, m, n] = b₂
            end
        end
    end

    @debug "Calculating D..."
    D₀₀ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₀₋₀ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₂₂ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₂₋₂ = OffsetArray(zeros(2N + 1, N, N), (-N):N, 1:N, 1:N)
    D₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1, N, N), (-N):N, 1:N, 1:N)

    Threads.@threads for n in 1:N
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

    @debug "Calculating g..."
    g₀₀ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₀₋₀ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₂₂ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₂₋₂ = OffsetArray(zeros(2N + 1), 0:(2N))
    g₀₂ = OffsetArray(zeros(ComplexF64, 2N + 1), 0:(2N))

    Threads.@threads for l in 0:(2N)
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

    @sync for i in 1:Threads.nthreads()
        StableTasks.@spawnat i + tid_offset wig_temp_free()
    end
    wig_table_free()

    return α₁, α₂, α₃, α₄, β₁, β₂
end

@testitem "Can calculate scattering matrix from axisymmetric T-Matrix" begin
    s = Spheroid(1.0, 2.0, complex(1.311))
    λ = 2π
    𝐓 = transition_matrix(s, λ)
    θs = collect(0:180)
    𝐅 = scattering_matrix(𝐓, λ, θs)

    @test size(𝐅) == (181, 6)
end

@testitem "Expansion coefficients can also be calculated by general function" begin
    s = Spheroid(1.0, 2.0, complex(1.311))
    λ = 2π
    𝐓 = transition_matrix(s, λ)
    T = TransitionMatrix{ComplexF64, size(𝐓, 2), typeof(𝐓)}(𝐓)

    α₁, α₂, α₃, α₄, β₁, β₂ = expansion_coefficients(𝐓, λ)
    α₁′, α₂′, α₃′, α₄′, β₁′, β₂′ = expansion_coefficients(T, λ)

    @test all(isapprox.(α₁, α₁′, atol = 1e-14))
    @test all(isapprox.(α₂, α₂′, atol = 1e-14))
    @test all(isapprox.(α₃, α₃′, atol = 1e-14))
    @test all(isapprox.(α₄, α₄′, atol = 1e-14))
    @test all(isapprox.(β₁, β₁′, atol = 1e-14))
    @test all(isapprox.(β₂, β₂′, atol = 1e-14))
end
