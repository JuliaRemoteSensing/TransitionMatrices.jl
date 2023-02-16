Base.convert(::Type{Float128}, x::ArbLike) = Float128(BigFloat(x))
Base.round(x::Arb, ::RoundingMode{:Up}) = ceil(BigFloat(x))
