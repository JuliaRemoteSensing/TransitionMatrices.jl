Base.convert(::Type{Float128}, x::ArbLike) = Float128(BigFloat(x))
