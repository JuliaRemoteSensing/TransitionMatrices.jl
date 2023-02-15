"""
Spherical Bessel function of the first kind.
"""
@inline function spherical_jn(::Type{T}, n::Integer, z) where {T}
    if T <: Arblib.ArbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_j!(Arblib.Arb(), Arblib.Arb(n + 1 // 2), Arblib.Arb(z))
    elseif T <: Arblib.AcbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_j!(Arblib.Acb(), Arblib.Acb(n + 1 // 2), Arblib.Acb(z))
    else
        return SpecialFunctions.sphericalbesselj(n, T(z))
    end
end

@inline spherical_jn(n::Integer, z::Real) = spherical_jn(Float64, n, z)
@inline spherical_jn(n::Integer, z) = spherical_jn(ComplexF64, n, z)

"""
First-order derivative of spherical Bessel function of the first kind.
"""
@inline function spherical_jn_deriv(::Type{T}, n::Integer, z) where {T}
    return spherical_jn(T, n - 1, z) - (n + 1) / z * spherical_jn(T, n, z)
end
@inline spherical_jn_deriv(n::Integer, z::Real) = spherical_jn_deriv(Float64, n, z)
@inline spherical_jn_deriv(n::Integer, z) = spherical_jn_deriv(ComplexF64, n, z)


"""
Spherical Bessel function of the second kind.
"""
function spherical_yn(::Type{T}, n::Integer, z) where {T}
    if T <: Arblib.ArbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_y!(Arblib.Arb(), Arblib.Arb(n + 1 // 2), Arblib.Arb(z))
    elseif T <: Arblib.AcbLike
        return √(Arblib.Arb(π) / (2z)) * Arblib.hypgeom_bessel_y!(Arblib.Acb(), Arblib.Acb(n + 1 // 2), Arblib.Acb(z))
    else
        return SpecialFunctions.sphericalbessely(n, T(z))
    end
end

@inline spherical_yn(n::Integer, z::Real) = spherical_yn(Float64, n, z)
@inline spherical_yn(n::Integer, z) = spherical_yn(ComplexF64, n, z)

"""
First-order derivative of spherical Bessel function of the second kind.
"""
@inline function spherical_yn_deriv(::Type{T}, n::Integer, z) where {T}
    return spherical_yn(T, n - 1, z) - (n + 1) / z * spherical_yn(T, n, z)
end
@inline spherical_yn_deriv(n::Integer, z::Real) = spherical_yn_deriv(Float64, n, z)
@inline spherical_yn_deriv(n::Integer, z) = spherical_yn_deriv(ComplexF64, n, z)

"""
Spherical Hankel function of the first kind.
"""
@inline spherical_hn1(::Type{T}, n::Integer, z) where {T} = spherical_jn(T, n, z) + 1im * spherical_yn(T, n, z)
@inline spherical_hn1(n::Integer, z::Real) = spherical_hn1(Float64, n, z)
@inline spherical_hn1(n::Integer, z) = spherical_hn1(ComplexF64, n, z)

"""
First-order derivative of spherical Hankel function of the first kind.
"""
@inline function spherical_hn1_deriv(::Type{T}, n::Integer, z) where {T}
    return spherical_jn_deriv(T, n, z) + 1im * spherical_yn_deriv(T, n, z)
end

@inline spherical_hn1_deriv(n::Integer, z::Real) = spherical_hn1_deriv(Float64, n, z)
@inline spherical_hn1_deriv(n::Integer, z) = spherical_hn1_deriv(ComplexF64, n, z)

"""
Spherical Hankel function of the second kind.
"""
@inline spherical_hn2(::Type{T}, n::Integer, z) where {T} = spherical_jn(T, n, z) - 1im * spherical_yn(T, n, z)

@inline spherical_hn2(n::Integer, z::Real) = spherical_hn2(Float64, n, z)
@inline spherical_hn2(n::Integer, z) = spherical_hn2(ComplexF64, n, z)

"""
First-order derivative of spherical Hankel function of the second kind.
"""
@inline function spherical_hn2_deriv(::Type{T}, n::Integer, z) where {T}
    return spherical_jn_deriv(T, n, z) - 1im * spherical_yn_deriv(T, n, z)
end

@inline spherical_hn2_deriv(n::Integer, z::Real) = spherical_hn2_deriv(Float64, n, z)
@inline spherical_hn2_deriv(n::Integer, z) = spherical_hn2_deriv(ComplexF64, n, z)

"""
Riccati Bessel function of the first kind.
"""
@inline ricatti_jn(::Type{T}, n::Integer, z) where {T} = z * spherical_jn(T, n, z)

@inline ricatti_jn(n::Integer, z::Real) = ricatti_jn(Float64, n, z)
@inline ricatti_jn(n::Integer, z) = ricatti_jn(ComplexF64, n, z)

"""
First-order derivative of Riccati Bessel function of the first kind.
"""
@inline function ricatti_jn_deriv(::Type{T}, n::Integer, z) where {T}
    jn = spherical_jn(T, n, z)
    jnd = spherical_jn_deriv(T, n, z)
    return z * jnd + jn
end

@inline ricatti_jn_deriv(n::Integer, z::Real) = ricatti_jn_deriv(Float64, n, z)
@inline ricatti_jn_deriv(n::Integer, z) = ricatti_jn_deriv(ComplexF64, n, z)

"""
Riccati Bessel function of the second kind.

> Note that in `miepy`, the author used `-z⋅y(z)` instead of `z⋅y(z)`
"""
@inline ricatti_yn(::Type{T}, n::Integer, z) where {T} = z * spherical_yn(T, n, z)

@inline ricatti_yn(n::Integer, z::Real) = ricatti_yn(Float64, n, z)
@inline ricatti_yn(n::Integer, z) = ricatti_yn(ComplexF64, n, z)

"""
First-order derivative of Riccati Bessel function of the second kind.
"""
@inline function ricatti_yn_deriv(::Type{T}, n::Integer, z) where {T}
    yn = spherical_yn(T, n, z)
    ynd = spherical_yn_deriv(T, n, z)
    return z * ynd + yn
end

@inline ricatti_yn_deriv(n::Integer, z::Real) = ricatti_yn_deriv(Float64, n, z)
@inline ricatti_yn_deriv(n::Integer, z) = ricatti_yn_deriv(ComplexF64, n, z)

"""
Riccati Hankel function of the first kind.
"""
@inline ricatti_hn1(::Type{T}, n::Integer, z) where {T} = ricatti_jn(T, n, z) + 1im * ricatti_yn(T, n, z)

@inline ricatti_hn1(n::Integer, z::Real) = ricatti_hn1(Float64, n, z)
@inline ricatti_hn1(n::Integer, z) = ricatti_hn1(ComplexF64, n, z)

"""
First-order derivative of Riccati Hankel function of the first kind.
"""
@inline function ricatti_hn1_deriv(::Type{T}, n::Integer, z) where {T}
    return ricatti_jn_deriv(T, n, z) + 1im * ricatti_yn_deriv(T, n, z)
end

@inline ricatti_hn1_deriv(n::Integer, z::Real) = ricatti_hn1_deriv(Float64, n, z)
@inline ricatti_hn1_deriv(n::Integer, z) = ricatti_hn1_deriv(ComplexF64, n, z)

"""
Riccati Hankel function of the second kind.
"""
@inline ricatti_hn2(::Type{T}, n::Integer, z) where {T} = ricatti_jn(T, n, z) - 1im * ricatti_yn(T, n, z)

@inline ricatti_hn2(n::Integer, z::Real) = ricatti_hn2(Float64, n, z)
@inline ricatti_hn2(n::Integer, z) = ricatti_hn2(ComplexF64, n, z)

"""
First-order derivative of Riccati Hankel function of the second kind.
"""
@inline function ricatti_hn2_deriv(::Type{T}, n::Integer, z) where {T}
    return ricatti_jn_deriv(T, n, z) - 1im * ricatti_yn_deriv(T, n, z)
end

@inline ricatti_hn2_deriv(n::Integer, z::Real) = ricatti_hn2_deriv(Float64, n, z)
@inline ricatti_hn2_deriv(n::Integer, z) = ricatti_hn2_deriv(ComplexF64, n, z)
