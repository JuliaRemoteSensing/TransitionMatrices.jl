@doc raw"""
```
factorial([T=Float64,], n)
```

Calculate factorials ``n!``. For ``n\le20``, the standard `Base.factorial` is used. For ``n>20``, `Arblib.gamma!` is used instead.
"""
function factorial(::Type{T}, n)::T where {T}
    if n <= 20 && isbits(T)
        return T(Base.factorial(n))
    end

    return T(Arblib.gamma!(Arb(), Arb(n + 1); prec = precision(T)))
end

@inline factorial(n) = factorial(Float64, n)
