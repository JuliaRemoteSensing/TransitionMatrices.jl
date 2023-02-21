const FACTORIAL = Dict()

"""
```
factorial([T=Float64,], n)
```

Calculate factorials with a global cache for every type.

!!! note
    This function is not thread-safe. To use it in a multi-threaded environment,
    you should pre-compute the maximum factorial you need for a specific type.

    For dynamic-precision types, e.g., `BigFloat` or `Arb`, after you change the 
    precision, you should call `clear_factorial_table!(T)` to clear the cache, and then
    re-compute the maximum factorial you need for the new precision.

    By default, the maximum factorial for the following types are pre-computed:
    - `Float64`: 150
    - `Double64`: 150
    - `Float128`: 300
    - `Arb`: 500
    - `BigFloat`: 500
"""
function factorial(::Type{T}, n)::T where {T}
    if n <= 20 && isbits(T)
        return T(Base.factorial(n))
    end

    if haskey(FACTORIAL, T)
        memo = FACTORIAL[T]
        n₀ = length(memo)
        if n₀ <= n
            sizehint!(memo, n + 1)
            for i in n₀:n
                push!(memo, memo[end] * i)
            end
        end
        return memo[n + 1]
    else
        FACTORIAL[T] = T[T(Base.factorial(i)) for i in 0:20]
        return factorial(T, n)
    end
end

@inline factorial(n) = factorial(Float64, n)

"""
```
clear_factorial_table!(T)
```

Clear the cached factorial table for a specific type. This is necessary after you 
change the precision of a dynamic-precision type, e.g., `Arb` or `BigFloat`.
"""
function clear_factorial_table!(::Type{T}) where {T}
    if haskey(FACTORIAL, T)
        pop!(FACTORIAL, T)
    end
end

@testitem "clear_factorial_table!" begin
    using TransitionMatrices: factorial, clear_factorial_table!

    clear_factorial_table!(BigFloat)
    setprecision(BigFloat, 300)
    @test precision(factorial(BigFloat, 500)) == 300

    setprecision(BigFloat, 256)
    clear_factorial_table!(BigFloat)
    @test precision(factorial(BigFloat, 500)) == 256
end
