Base.convert(::Type{Float128}, x::ArbLike) = Float128(BigFloat(x))
Quadmath.Float128(x::ArbLike) = Float128(BigFloat(x))
Base.round(x::Arb, ::RoundingMode{:Up}) = ceil(BigFloat(x))
Base.abs2(x::AcbLike) = abs2(real(x)) + abs2(imag(x))
Base.complex(::Type{Arb}) = Acb

# Workaround for DoubleFloats v1.9.x: `cbrt` (and thus its alias `∛`) on a
# `Double64` returns the raw `(hi, lo)` component tuple instead of a `Double64`
# (the upstream `cbrt_db_db` forgets to wrap its result; `sqrt` is unaffected).
# This breaks every `∛` call reached with a `Double64` argument, e.g.
# `volume_equivalent_radius` and the Riccati-Bessel term estimators.
#
# The fix is behavior-gated: the guard below calls the upstream `cbrt` before
# our method exists, so it installs the workaround only while upstream is still
# broken. Once a fixed DoubleFloats returns a `Double64` here, this block is
# skipped at load time and the upstream method is used unchanged — no sticky
# shadowing, no dead code. (A package upgrade recompiles this module, so the
# guard is re-evaluated against the new DoubleFloats version.)
if !(cbrt(Double64(8.0)) isa Double64)
    # Specialized to `Double64` (more specific than upstream's parametric
    # method); recovers full precision with two Newton steps from a Float64
    # seed. `cbrt(Float64(a))` dispatches to Base, so there is no recursion.
    function Base.cbrt(x::Double64)
        iszero(x) && return x
        a = abs(x)
        y = Double64(cbrt(Float64(a)))
        y = (2y + a / (y * y)) / 3
        y = (2y + a / (y * y)) / 3
        return x < 0 ? -y : y
    end
end

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
