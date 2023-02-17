Base.convert(::Type{Float128}, x::ArbLike) = Float128(BigFloat(x))
Base.round(x::Arb, ::RoundingMode{:Up}) = ceil(BigFloat(x))

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
