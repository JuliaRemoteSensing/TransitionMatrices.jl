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

# Compat for enumerate
Base.firstindex(::Base.Iterators.Enumerate{OrderDegreeIterator}) = 1
Base.lastindex(iter::Base.Iterators.Enumerate{OrderDegreeIterator}) = length(iter.itr)
function Base.getindex(iter::Base.Iterators.Enumerate{OrderDegreeIterator}, idx)
    (idx, iter.itr[idx])
end
