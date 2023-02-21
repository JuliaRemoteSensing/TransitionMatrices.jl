struct RandomOrientationTransitionMatrix{CT,N,V<:AbstractArray{CT, 3}} <: AbstractTransitionMatrix{CT,N}
    𝐓::V
end

@doc raw"""
```
RandomOrientationTransitionMatrix(𝐓::AbstractTransitionMatrix{CT, N}) where {CT, N}
```

Calculate the random-orientation T-Matrix according to Eq. (5.96) in Mishchenko et al. (2002).

```math
\left\langle T_{m n m^{\prime} n^{\prime}}^{k l}(L)\right\rangle=\frac{\delta_{n n^{\prime}} \delta_{m m^{\prime}}}{2 n+1} \sum_{m_1=-n}^n T_{m_1 n m_1 n}^{k l}(P), \quad k, l=1,2
```
"""
function RandomOrientationTransitionMatrix(𝐓::AbstractTransitionMatrix{CT, N}) where {CT, N}
    T̄ = OffsetArray(zeros(CT, N, 2, 2), 1:N, 1:2, 1:2)
    @inbounds for p′ in 1:2, p in 1:2, n in 1:N
        T̄[n, p, p′] = sum(𝐓[m₁, n, m₁, n, p, p′] for m₁ in -n:n) / (2n+1)
    end
    return RandomOrientationTransitionMatrix{CT,N,typeof(T̄)}(T̄)
end

Base.@propagate_inbounds function Base.getindex(oa::RandomOrientationTransitionMatrix{CT,
                                                                                        N, V
                                                                                        },
                                                m::Integer, n::Integer, m′::Integer,
                                                n′::Integer, p::Integer,
                                                p′::Integer) where {CT, N, V}
    if m != m′ || n != n′ || abs(m) > n
        zero(CT)
    else
        oa.𝐓[n, p, p′]
    end
end

rotate(ro::RandomOrientationTransitionMatrix{CT, N, V}, ::Rotation{3}) where {CT, N, V} = ro

orientation_average(ro::RandomOrientationTransitionMatrix, _pₒ; _kwargs...) = ro

@testitem "Analytical orientation average can be approximated by numerical method" begin
    using TransitionMatrices

    s = Spheroid{Float64, ComplexF64}(2, 4, 1.311);
    𝐓 = transition_matrix(s, 2π, 5, 20);

    T̄ = orientation_average(𝐓, (α, β, γ) -> 1 / 8π^2; Nα=40, Nβ=40, Nγ=1)
    T̄′ = RandomOrientationTransitionMatrix(𝐓)
    @test all(isapprox.(T̄, T̄′; atol=1e-14))
end
