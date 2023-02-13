"""
Iterator for the order-degree pairs of the given maximum order `nₘₐₓ`.
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
function Base.getindex(iter::OrderDegreeIterator, idx)
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

    @test iterate(OrderDegreeIterator(3)) == ((1, -1), (1, -1))
    @test iterate(OrderDegreeIterator(3), (1, 1)) == ((2, -2), (2, -2))
    @test collect(OrderDegreeIterator(2)) ==
          [(1, -1), (1, 0), (1, 1), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
    @test size(OrderDegreeIterator(100)) == (10200,)
    @test !Base.isdone(OrderDegreeIterator(2), (1, 1))
    @test Base.isdone(OrderDegreeIterator(2), (2, 2))
end

@doc raw"""
A general T-Matrix ``T_{m n m^{\prime} n^{\prime}}^{k l}`` stored in a 6-dimensional array, in the order ``(m, n, m^{\prime}, n^{\prime}, k, l)``.
"""
abstract type AbstractTransitionMatrix{CT, N} <: AbstractArray{CT, 6} end

"""
Get the maximum order of a T-Matrix.
"""
order(::AbstractTransitionMatrix{CT, N}) where {CT, N} = N

Base.size(::AbstractTransitionMatrix{CT, N}) where {CT, N} = (2N + 1, N, 2N + 1, N, 2, 2)
function Base.axes(::AbstractTransitionMatrix{CT, N}) where {CT, N}
    ((-N):N, 1:N, (-N):N, 1:N, 1:2, 1:2)
end

"""
Concrete type for a general T-Matrix.
"""
struct TransitionMatrix{CT, N, V <: AbstractArray{CT}} <: AbstractTransitionMatrix{CT, N}
    container::V
end

Base.getindex(tm::TransitionMatrix{CT, N}, idx) where {CT, N} = getindex(tm.container, idx)
function Base.getindex(tm::TransitionMatrix{CT, N}, idxs...) where {CT, N}
    getindex(tm.container, idxs...)
end

@doc raw"""
`rotate(𝐓::AbstractTransitionMatrix{CT, N}, rot::Rotation{3})`

Rotate the given T-Matrix `𝐓` by the Euler angle `rot` and generate a new T-Matrix.

- For a general T-Matrix, Eq. (5.29) in Mishchenko et al. (2002) is used as a fallback. A `TransitionMatrix` will be returned, which is the most general yet concrete type.

```math
T_{m n m^{\prime} n^{\prime}}^{p p′}(L ; \alpha, \beta, \gamma)=\sum_{m_1=-n}^n \sum_{m_2=-n^{\prime}}^{n^{\prime}} D_{m m_1}^n(\alpha, \beta, \gamma) T_{m_1 n m_2 n^{\prime}}^{p p′}(P) D_{m_2 m^{\prime}}^{n^{\prime}}(-\gamma,-\beta,-\alpha)\quad p,p′=1,2
```

- For a `MieTransitionMatrix`, the underlying Mie coefficients are copied and a new `MieTransitionMatrix` will be returned.
"""
function rotate(𝐓::AbstractTransitionMatrix{CT, N}, rot::Rotation{3}) where {CT, N}
    # Get the Euler angle in Z-Y-Z order.
    zyz = RotZYZ(rot)
    α, β, γ = zyz.theta1, zyz.theta2, zyz.theta3

    # Calculate the wigner-d functions that will be used.
    d = OffsetArray(zeros(CT, 2N + 1, 2N + 1, N + 1), (-N):N, (-N):N, 0:N)
    for m in (-N):N
        for m′ in (-N):N
            sₘᵢₙ = max(abs(m), abs(m′))
            wigner_d_recursion!(view(d, m, m′, sₘᵢₙ:N), m, m′, N, β)
        end
    end

    # Calculate the coefficients used for wigner-D functions
    coeff = OffsetArray([cis(-(m * α + m′ * γ)) for m in (-N):N, m′ in (-N):N], (-N):N,
                        (-N):N)

    # Calculate the rotated T-Matrix
    𝐓′ = similar(𝐓)
    fill!(𝐓′, 0)

    # Enable multi-threading
    Threads.@threads for (n′, m′) in OrderDegreeIterator(N)
        for p in 1:2, p′ in 1:2
            for (n, m) in OrderDegreeIterator(N)
                for m₂ in (-n′):n′, m₁ in (-n):n
                    sign = iseven(m′ + m₂) ? 1 : -1
                    𝐓′[m, n, m′, n′, p, p′] += coeff[m, m₁] * d[m, m₁, n] *
                                               conj(coeff[m′, m₂]) * d[m₂, m′, n′] * sign *
                                               𝐓[m₁, n, m₂, n′, p, p′]
                end
            end
        end
    end

    TransitionMatrix{CT, N, typeof(𝐓′)}(𝐓′)
end
