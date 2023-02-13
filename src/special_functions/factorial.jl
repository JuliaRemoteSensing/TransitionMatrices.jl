const FACTORIAL = Dict()

"""
Calculate factorials using the Gamma function.
"""
function factorial(::Type{T}, n)::T where {T}
    if n <= 20
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
