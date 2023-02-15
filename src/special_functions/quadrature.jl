gausslegendre(::Type{Float64}, n::Integer) = FGQ.gausslegendre(n)

function gausslegendre(T, n::Integer)
    x = ArbRefVector(n)
    w = ArbRefVector(n)
    gausslegendre!(T, n, x, w)
    x = convert.(T, x)
    w = convert.(T, w)
    return x, w
end

function gausslegendre(T::Type{<:Arblib.ArbLike}, n::Integer)
    x = ArbRefVector(n)
    w = ArbRefVector(n)
    gausslegendre!(T, n, x, w)
    return x, w
end

function gausslegendre!(::Type{<:Arblib.ArbLike}, n::Integer, x::Arblib.ArbVectorLike,
                        w::Arblib.ArbVectorLike)
    for i in 1:(n ÷ 2)
        Arblib.hypgeom_legendre_p_ui_root!(x[n + 1 - i], w[n + 1 - i], UInt64(n),
                                           UInt64(i - 1))
        x[i] = -x[n + 1 - i]
        w[i] = w[n + 1 - i]
    end

    if n % 2 == 1
        Arblib.hypgeom_legendre_p_ui_root!(x[n ÷ 2 + 1], w[n ÷ 2 + 1], UInt64(n),
                                           UInt64(n ÷ 2))
    end
end