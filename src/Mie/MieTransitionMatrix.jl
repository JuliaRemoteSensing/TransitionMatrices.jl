@doc raw"""
According to Eq. (5.42) -- Eq. (5.44) in Mishchenko et al. (2002), the T-Matrix for a Mie particle can be written as:

```math
\begin{array}{l}
T_{m n m^{\prime} n^{\prime}}^{12}(P) \equiv 0, \quad T_{m n m^{\prime} n^{\prime}}^{21}(P) \equiv 0, \\
T_{m n m^{\prime} n^{\prime}}^{11}(P)=-\delta_{m m^{\prime}} \delta_{n n^{\prime}} b_n, \\
T_{m n m^{\prime} n^{\prime}}^{22}(P)=-\delta_{m m^{\prime}} \delta_{n n^{\prime}} a_n .
\end{array}
```

```
MieTransitionMatrix{CT, N}(x::Real, m::Number)
``` 

Generate the T-Matrix from the Mie coefficients of a homogeneous sphere.

```
MieTransitionMatrix{CT, N}(x_core::Real, x_mantle::Real, m_core::Number, m_mantle::Number)
```

Generate the T-Matrix from the Mie coefficients of a coated sphere.

This struct provides the T-Matrix API for a Mie particle.
"""
struct MieTransitionMatrix{CT, N, V <: AbstractVector{CT}} <:
       AbstractTransitionMatrix{CT, N}
    a::V
    b::V
end

function MieTransitionMatrix{CT, N}(x::Real, m::Number) where {CT, N}
    T = real(CT)
    a, b = bhmie(T, x, m; nₘₐₓ = N)
    MieTransitionMatrix{CT, N, Vector{CT}}(a, b)
end

function MieTransitionMatrix{CT, N}(x_core::Real, x_mantle::Real, m_core::Number,
        m_mantle::Number) where {CT, N}
    T = real(CT)
    a, b = bhcoat(T, x_core, x_mantle, m_core, m_mantle; nₘₐₓ = N)
    MieTransitionMatrix{CT, N, Vector{CT}}(a, b)
end

Base.@propagate_inbounds function Base.getindex(mie::MieTransitionMatrix{CT, N, V},
        m::Integer, n::Integer, m′::Integer,
        n′::Integer, p::Integer,
        p′::Integer) where {CT, N, V}
    if m != m′ || n != n′ || p != p′ || abs(m) > n
        zero(CT)
    else
        p == 1 ? -mie.b[n] : -mie.a[n]
    end
end

"""
```
rotate(mie::MieTransitionMatrix, ::Rotation{3})
```

The T-Matrix of a Mie scatterer is invariant under rotation. Therefore, the original T-Matrix will be returned.
"""
rotate(mie::MieTransitionMatrix, ::Rotation{3}) = mie

"""
```
orientation_average(mie::MieTransitionMatrix, _pₒ; _kwargs...)
```

The T-Matrix of a Mie scatterer is invariant under rotation. Therefore, the original T-Matrix will be returned.
"""
orientation_average(mie::MieTransitionMatrix, _pₒ; _kwargs...) = mie

function scattering_cross_section(mie::MieTransitionMatrix{CT, N, V},
        λ = 2π) where {CT, N, V}
    Cˢᶜᵃ = zero(real(CT))

    @inbounds for n in 1:N
        Cˢᶜᵃ += (2n + 1) * (abs2(mie.a[n]) + abs2(mie.b[n]))
    end

    return Cˢᶜᵃ * λ^2 / 2π
end

function extinction_cross_section(mie::MieTransitionMatrix{CT, N, V},
        λ = 2π) where {CT, N, V}
    Cᵉˣᵗ = zero(real(CT))

    @inbounds for n in 1:N
        Cᵉˣᵗ += (2n + 1) * real(mie.a[n] + mie.b[n])
    end

    return Cᵉˣᵗ * λ^2 / 2π
end

function amplitude_matrix(mie::MieTransitionMatrix{CT, N, V}, ϑᵢ, φᵢ, ϑₛ, φₛ;
        λ = 2π) where {CT, N, V}
    T = real(CT)
    k₁ = 2π / λ
    𝐒₁₁, 𝐒₁₂, 𝐒₂₁, 𝐒₂₂ = zero(CT), zero(CT), zero(CT), zero(CT)

    ϑᵢ = T(ϑᵢ)
    φᵢ = T(φᵢ)
    ϑₛ = T(ϑₛ)
    φₛ = T(φₛ)
    Δφ = φₛ - φᵢ

    πᵢ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    τᵢ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    πₛ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)
    τₛ = OffsetArray(zeros(T, 2N + 1, N + 1), (-N):N, 0:N)

    for m in 0:N
        wigner_d_recursion!(view(πᵢ, m, m:N), 0, m, N, ϑᵢ;
            deriv = view(τᵢ, m, m:N))
        wigner_d_recursion!(view(πₛ, m, m:N), 0, m, N, ϑₛ;
            deriv = view(τₛ, m, m:N))
    end

    for n in 1:N
        for m in 0:n
            πᵢ[m, n] = pi_func(T, m, n, ϑᵢ; d = πᵢ[m, n])
            πₛ[m, n] = pi_func(T, m, n, ϑₛ; d = πₛ[m, n])
            if m > 0
                π_sign = iseven(m + 1) ? one(T) : -one(T)
                τ_sign = iseven(m) ? one(T) : -one(T)
                πᵢ[-m, n] = π_sign * πᵢ[m, n]
                πₛ[-m, n] = π_sign * πₛ[m, n]
                τᵢ[-m, n] = τ_sign * τᵢ[m, n]
                τₛ[-m, n] = τ_sign * τₛ[m, n]
            end
        end
    end

    @inbounds for n in 1:N
        T₁₁ = -mie.b[n]
        T₂₂ = -mie.a[n]
        αₙ = CT(-1im) * T(2n + 1) / T(n * (n + 1))

        for m in (-n):n
            expiφ = cis(m * Δφ)
            πᵢₘₙ = πᵢ[m, n]
            τᵢₘₙ = τᵢ[m, n]
            πₛₘₙ = πₛ[m, n]
            τₛₘₙ = τₛ[m, n]

            𝐒₁₁ += αₙ * (T₁₁ * πₛₘₙ * πᵢₘₙ + T₂₂ * τₛₘₙ * τᵢₘₙ) * expiφ
            𝐒₁₂ += αₙ * (T₁₁ * πₛₘₙ * τᵢₘₙ + T₂₂ * τₛₘₙ * πᵢₘₙ) * expiφ
            𝐒₂₁ += αₙ * (T₁₁ * τₛₘₙ * πᵢₘₙ + T₂₂ * πₛₘₙ * τᵢₘₙ) * expiφ
            𝐒₂₂ += αₙ * (T₁₁ * τₛₘₙ * τᵢₘₙ + T₂₂ * πₛₘₙ * πᵢₘₙ) * expiφ
        end
    end

    return (@SMatrix [𝐒₁₁ 𝐒₁₂/1im; 𝐒₂₁*1im 𝐒₂₂]) ./ k₁
end

@testitem "MieTransitionMatrix" begin
    using TransitionMatrices: MieTransitionMatrix, RotZYZ, TransitionMatrix, rotate

    @testset "of a homogeneous sphere remains the same under rotations" begin
        𝐓 = MieTransitionMatrix{ComplexF64, 5}(1.0, 1.311)
        @test all(isapprox.(𝐓,
            rotate(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓),
                RotZYZ(0.2, 0.3, 0.4)); atol = eps(Float64)))
        @test all(isapprox.(𝐓,
            rotate(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓),
                RotZYZ(0.8, 0.0, -1.0)); atol = eps(Float64)))
    end

    @testset "of a coated sphere remains the same under rotations" begin
        𝐓 = MieTransitionMatrix{ComplexF64, 5}(0.4, 0.8, 1.0, 1.311)
        @test all(isapprox.(𝐓,
            rotate(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓),
                RotZYZ(0.2, 0.3, 0.4)); atol = eps(Float64)))
        @test all(isapprox.(𝐓,
            rotate(TransitionMatrix{ComplexF64, 5, typeof(𝐓)}(𝐓),
                RotZYZ(0.8, 0.0, -1.0)); atol = eps(Float64)))
    end
end

@testitem "MieTransitionMatrix uses specialized formulas consistent with fallback" begin
    using TransitionMatrices: MieTransitionMatrix, TransitionMatrix, amplitude_matrix,
                              extinction_cross_section, scattering_cross_section

    𝐓 = MieTransitionMatrix{ComplexF64, 6}(1.7, 1.311 + 0.01im)
    𝐓_fallback = TransitionMatrix{ComplexF64, 6, typeof(𝐓)}(𝐓)
    angles = (0.2, 0.3, 1.2, 0.5)
    λ = 2π

    @test scattering_cross_section(𝐓, λ) ≈ scattering_cross_section(𝐓_fallback, λ)
    @test extinction_cross_section(𝐓, λ) ≈ extinction_cross_section(𝐓_fallback, λ)
    @test all(isapprox.(amplitude_matrix(𝐓, angles...; λ = λ),
        amplitude_matrix(𝐓_fallback, angles...; λ = λ)))
end
