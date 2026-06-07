module TransitionMatricesQuadmathExt

# Float128 (Quadmath) support, loaded on demand. Quadmath is a *weak* dependency: the
# core package never needs Float128 (it is a purely user-opt-in numeric type), so the
# Float128 ⇄ Arb glue lives here and activates only when the user has `using Quadmath`.
# (Arblib is a hard dependency, so it is available to this extension.)

using TransitionMatrices
using Quadmath: Quadmath, Float128, ComplexF128
using Arblib: Arblib, ArbLike

Base.convert(::Type{Float128}, x::ArbLike) = Float128(BigFloat(x))
Quadmath.Float128(x::ArbLike) = Float128(BigFloat(x))
Base.precision(::Type{ComplexF128}) = 113

# Quadmath implements `precision(::Type{Float128})` but not the *instance* method
# `precision(::Float128)` (it routes to a missing `_precision_with_base_2`). Define it so
# generic `precision(x)` works on Float128 values (Float128 has a 113-bit significand).
Base.precision(::Float128) = 113

function Arblib.set!(arb::Arblib.ArbLike, val::Float128)
    Arblib.set!(arb, BigFloat(val))
end

end
