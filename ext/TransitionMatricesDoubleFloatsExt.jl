module TransitionMatricesDoubleFloatsExt

# Double64 (DoubleFloats) support, loaded on demand. DoubleFloats is a *weak* dependency:
# the core package never needs Double64 (it is a purely user-opt-in numeric type), so the
# Double64-specific shims live here and activate only when the user has `using DoubleFloats`.
# (Arblib is a hard dependency, so it is available to this extension.)

using TransitionMatrices
using DoubleFloats: Double64
using Arblib: Arblib

Base.precision(::Type{Complex{Double64}}) = 106

function Arblib.set!(arb::Arblib.ArbLike, val::Double64)
    Arblib.set!(arb, BigFloat(val))
end

# Workaround for DoubleFloats v1.9.x: `cbrt` (and thus its alias `∛`) on a `Double64`
# returns the raw `(hi, lo)` component tuple instead of a `Double64` (the upstream
# `cbrt_db_db` forgets to wrap its result; `sqrt` is unaffected). This breaks every `∛`
# call reached with a `Double64` argument, e.g. `volume_equivalent_radius` and the
# Riccati-Bessel term estimators.
#
# The fix is behavior-gated: the guard below calls the upstream `cbrt` before our method
# exists, so it installs the workaround only while upstream is still broken. Once a fixed
# DoubleFloats returns a `Double64` here, this block is skipped at load time and the
# upstream method is used unchanged — no sticky shadowing, no dead code.
if !(cbrt(Double64(8.0)) isa Double64)
    # Specialized to `Double64` (more specific than upstream's parametric method);
    # recovers full precision with two Newton steps from a Float64 seed. `cbrt(Float64(a))`
    # dispatches to Base, so there is no recursion.
    function Base.cbrt(x::Double64)
        iszero(x) && return x
        a = abs(x)
        y = Double64(cbrt(Float64(a)))
        y = (2y + a / (y * y)) / 3
        y = (2y + a / (y * y)) / 3
        return x < 0 ? -y : y
    end
end

end
