Base.convert(::Type{Float128}, x::ArbLike) = Float128(BigFloat(x))
Quadmath.Float128(x::ArbLike) = Float128(BigFloat(x))
Base.round(x::Arb, ::RoundingMode{:Up}) = ceil(BigFloat(x))
Base.abs2(x::AcbLike) = abs2(real(x)) + abs2(imag(x))

function Base.inv(x::Matrix{Arb})
    a = ArbMatrix(x)
    Arblib.inv!(a, a)
    a
end

function Base.inv(x::Matrix{Acb})
    a = AcbMatrix(x)
    Arblib.inv!(a, a)
    a
end

Base.precision(::Type{ComplexF32}) = 24
Base.precision(::Type{ComplexF64}) = 53
Base.precision(::Type{Complex{Double64}}) = 106
Base.precision(::Type{ComplexF128}) = 113
Base.precision(::Type{Complex{Arb}}) = precision(Arb)

function Base.precision(::Type{
                               Complex{ForwardDiff.Dual{ForwardDiff.Tag{F, T},
                                                        T, N}}}) where {F, T, N}
    precision(T)
end

function Arblib.set!(arb::Arblib.ArbLike, dual::ForwardDiff.Dual)
    Arblib.set!(arb, dual.value)
end

function Arblib.set!(arb::Arblib.ArbLike, val::Float128)
    Arblib.set!(arb, BigFloat(val))
end

function Arblib.set!(arb::Arblib.ArbLike, val::Double64)
    Arblib.set!(arb, BigFloat(val))
end

Base.Int64(x::ArbLike) = round(Int64, BigFloat(x))

function Base.Float16(x::T) where {T <: ForwardDiff.Dual}
    Float16(ForwardDiff.value(x))
end

Base.Float32(x::ArbLike) = Float32(BigFloat(x))
function Base.Float32(x::T) where {T <: ForwardDiff.Dual}
    Float32(ForwardDiff.value(x))
end

function Base.BigFloat(x::T) where {T <: ForwardDiff.Dual}
    BigFloat(ForwardDiff.value(x))
end

function Base.convert(::Type{Complex{ForwardDiff.Dual{ForwardDiff.Tag{F, T}, T, N}}},
                      x::AcbLike) where {F, T, N}
    re = ForwardDiff.Dual{ForwardDiff.Tag{F, T}, T, N}(real(x))
    im = ForwardDiff.Dual{ForwardDiff.Tag{F, T}, T, N}(imag(x))
    return Complex{ForwardDiff.Dual{ForwardDiff.Tag{F, T}, T, N}}(re, im)
end

function Base.convert(::Type{Complex{Arb}}, x::AcbLike)
    re = Arb(real(x))
    im = Arb(imag(x))
    return Complex{Arb}(re, im)
end
