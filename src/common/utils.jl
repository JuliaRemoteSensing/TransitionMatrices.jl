"""
Iterator for the order-degree pairs of the given maximum order `nₘₐₓ`.

Example of `nₘₐₓ=2`:

```jldoctest
julia> collect(OrderDegreeIterator(2))
8-element Vector{Tuple{Int64, Int64}}:
 (1, -1)
 (1, 0)
 (1, 1)
 (2, -2)
 (2, -1)
 (2, 0)
 (2, 1)
 (2, 2)
```
"""
struct OrderDegreeIterator
    nₘₐₓ::Int
end

Base.iterate(::OrderDegreeIterator) = ((1, -1), (1, -1))
function Base.iterate(iter::OrderDegreeIterator, (n, m))
    if m == n
        n == iter.nₘₐₓ ? nothing : ((n + 1, -n - 1), (n + 1, -n - 1))
    else
        ((n, m + 1), (n, m + 1))
    end
end

Base.firstindex(iter::OrderDegreeIterator) = 1
Base.lastindex(iter::OrderDegreeIterator) = length(iter)
function Base.getindex(::OrderDegreeIterator, idx)
    n = floor(Int, √idx)
    m = idx - n^2 - n
    (n, m)
end
Base.length(iter::OrderDegreeIterator) = iter.nₘₐₓ * (iter.nₘₐₓ + 2)
Base.size(x::OrderDegreeIterator) = (length(x),)
Base.eltype(::OrderDegreeIterator) = Tuple{Int, Int}
Base.isdone(iter::OrderDegreeIterator, state) = state >= (iter.nₘₐₓ, iter.nₘₐₓ)

@testitem "OrderDegreeIterator" begin
    using TransitionMatrices: OrderDegreeIterator

    @test eltype(OrderDegreeIterator(3)) == Tuple{Int, Int}
    @test iterate(OrderDegreeIterator(3)) == ((1, -1), (1, -1))
    @test iterate(OrderDegreeIterator(3), (1, 1)) == ((2, -2), (2, -2))
    @test collect(OrderDegreeIterator(2)) ==
          [(1, -1), (1, 0), (1, 1), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
    @test length(OrderDegreeIterator(2)) == 8
    @test size(OrderDegreeIterator(100)) == (10200,)
    @test !Base.isdone(OrderDegreeIterator(2), (1, 1))
    @test Base.isdone(OrderDegreeIterator(2), (2, 2))
end

"""
```
phase_matrix(𝐒::AbstractMatrix)
```

Calculate the phase matrix `𝐙` from the amplitude matrix `𝐒`, according to Eq. (2.106) -- Eq. (2.121) in Mishchenko et al. (2002).
"""
function phase_matrix(𝐒::AbstractMatrix)
    𝐙₁₁ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' + 𝐒[1, 2] * 𝐒[1, 2]' + 𝐒[2, 1] * 𝐒[2, 1]' +
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₁₂ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[1, 2]' + 𝐒[2, 1] * 𝐒[2, 1]' -
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₁₃ = -𝐒[1, 1] * 𝐒[1, 2]' - 𝐒[2, 2] * 𝐒[2, 1]'
    𝐙₁₄ = 1.0im * (𝐒[1, 1] * 𝐒[1, 2]' - 𝐒[2, 2] * 𝐒[2, 1]')
    𝐙₂₁ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' + 𝐒[1, 2] * 𝐒[1, 2]' - 𝐒[2, 1] * 𝐒[2, 1]' -
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₂₂ = 0.5 * (𝐒[1, 1] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[1, 2]' - 𝐒[2, 1] * 𝐒[2, 1]' +
           𝐒[2, 2] * 𝐒[2, 2]')
    𝐙₂₃ = -𝐒[1, 1] * 𝐒[1, 2]' + 𝐒[2, 2] * 𝐒[2, 1]'
    𝐙₂₄ = 1.0im * (𝐒[1, 1] * 𝐒[1, 2]' + 𝐒[2, 2] * 𝐒[2, 1]')
    𝐙₃₁ = -𝐒[1, 1] * 𝐒[2, 1]' - 𝐒[2, 2] * 𝐒[1, 2]'
    𝐙₃₂ = -𝐒[1, 1] * 𝐒[2, 1]' + 𝐒[2, 2] * 𝐒[1, 2]'
    𝐙₃₃ = 𝐒[1, 1] * 𝐒[2, 2]' + 𝐒[1, 2] * 𝐒[2, 1]'
    𝐙₃₄ = -1.0im * (𝐒[1, 1] * 𝐒[2, 2]' + 𝐒[2, 1] * 𝐒[1, 2]')
    𝐙₄₁ = 1.0im * (𝐒[2, 1] * 𝐒[1, 1]' + 𝐒[2, 2] * 𝐒[1, 2]')
    𝐙₄₂ = 1.0im * (𝐒[2, 1] * 𝐒[1, 1]' - 𝐒[2, 2] * 𝐒[1, 2]')
    𝐙₄₃ = -1.0im * (𝐒[2, 2] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[2, 1]')
    𝐙₄₄ = 𝐒[2, 2] * 𝐒[1, 1]' - 𝐒[1, 2] * 𝐒[2, 1]'

    𝐙 = @SMatrix [𝐙₁₁ 𝐙₁₂ 𝐙₁₃ 𝐙₁₄; 𝐙₂₁ 𝐙₂₂ 𝐙₂₃ 𝐙₂₄; 𝐙₃₁ 𝐙₃₂ 𝐙₃₃ 𝐙₃₄; 𝐙₄₁ 𝐙₄₂ 𝐙₄₃ 𝐙₄₄]
    return real.(𝐙)
end

@testitem "Can calculate phase matrix from amplitude scattering matrix" begin
    using TransitionMatrices

    @test all(phase_matrix([1+2im 2+3im; 0.2-0.5im 0.5-0.2im]) .≈
              [9.29 -4.0 -8.2 -0.79
               8.71 -4.0 -7.8 -1.21
               0.4 1.2 -1.0 -0.4
               2.8 -1.0 -2.8 1.2])
end

"""
```
albedo(𝐓::AbstractTransitionMatrix)
```

Calculate the single scattering albedo from the given T-Matrix.
"""
function albedo(𝐓::AbstractTransitionMatrix)
    scattering_cross_section(𝐓) / extinction_cross_section(𝐓)
end

@testitem "albedo of" begin
    using TransitionMatrices: Spheroid, transition_matrix, albedo

    @testset "non-absorbing scatterers should be equal to 1.0" begin
        s = Spheroid(1.5, 1.0, complex(1.311))
        𝐓 = transition_matrix(s, 2π, 5, 40)
        @test albedo(𝐓) ≈ 1.0
    end

    @testset "absorbing scatterers should be less than 1.0" begin
        s = Spheroid(1.5, 1.0, 1.5 + 0.01im)
        𝐓 = transition_matrix(s, 2π, 5, 40)
        @test albedo(𝐓) < 1.0
    end
end

"""
```
absorption_cross_section(𝐓::AbstractTransitionMatrix, λ=2π)
```

Calculate the absorption cross section from the given T-Matrix.
"""
function absorption_cross_section(𝐓::AbstractTransitionMatrix, λ = 2π)
    extinction_cross_section(𝐓, λ) - scattering_cross_section(𝐓, λ)
end

@testitem "absorption cross section of" begin
    using TransitionMatrices: Spheroid, transition_matrix, absorption_cross_section

    @testset "non-absorbing scatterers should be approximate to 0.0" begin
        s = Spheroid(1.5, 1.0, complex(1.311))
        𝐓 = transition_matrix(s, 2π, 5, 40)

        # The result may be slightly negative due to numerical errors
        @test abs(absorption_cross_section(𝐓)) < 1e-8
    end

    @testset "absorbing scatterers should be less than 1.0" begin
        s = Spheroid(1.5, 1.0, 1.5 + 0.01im)
        𝐓 = transition_matrix(s, 2π, 5, 40)
        @test absorption_cross_section(𝐓) > 0.0
    end
end
