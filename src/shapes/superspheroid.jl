@doc raw"""
A super-spheroid scatterer (Bi/Lin dust model, D4h symmetry).

Defined by the implicit surface

```math
\left|\frac{x}{a}\right|^{2/n} + \left|\frac{y}{a}\right|^{2/n} + \left|\frac{z}{c}\right|^{2/n} = 1
```

where the cross-sections at constant ``z`` are super-circles (rounded squares for
``n \neq 1``), giving **4-fold** (D4h) rotational symmetry — NOT axial symmetry.
Use `SuperSpheroidRevolved` for the true surface of revolution.

Conventions (Lin & Bi 2018 JGR, doi:10.1029/2018JD029464):
- ``n = 1``: ordinary spheroid ``(a = b)``.
- ``n = 2``: regular octahedron.
- ``1 < n < 2``: convex rounded shape.
- ``n > 2``: concave (dust best-fit region).

Attributes:

- `a`: the equatorial semi-axis (``x`` and ``y`` directions).
- `c`: the polar semi-axis (``z`` direction).
- `n`: the roundness exponent.
- `m`: the relative complex refractive index.
"""
struct SuperSpheroid{T, CT} <: AbstractNFoldShape{4, T, CT}
    a::T
    c::T
    n::T
    m::CT
end

SuperSpheroid(a, c, n, m) = SuperSpheroid{typeof(a), typeof(m)}(a, c, n, m)

# Beta function B(x, y) = Γ(x)Γ(y)/Γ(x+y) for the superquadric volume formulas.
# Computed in Arb (lgamma) at the inputs' own precision and converted straight back to
# `T` — the same `T(Arblib...(Arb(...); prec))` idiom the package uses elsewhere (see
# `special_functions/factorial.jl`). No BigFloat intermediate: that would both be a
# needless detour and, because `BigFloat(::Arb)` uses the *global* precision, silently
# truncate high-precision `T` (Arb, large BigFloat). `precision(x)` reads the instance's
# precision, so Float64/Double64/Float128/BigFloat/Arb shapes each compute at their own.
# Working precision (bits) from the value's own precision, so high-precision shapes
# (Double64, BigFloat, Arb instances) are not truncated. `precision(::Float128)` (the
# instance method) is unimplemented in the current Quadmath/Julia (it routes to a
# missing `_precision_with_base_2`), so Float128 is pinned to its 113-bit mantissa via
# the working type method `precision(Float128)`.
_beta_prec(::Float128) = precision(Float128)
_beta_prec(x::Real) = precision(x)

function _beta_lgamma(x::T, y::T) where {T <: Real}
    prec = max(_beta_prec(x), _beta_prec(y))
    lΓx = Arblib.lgamma!(Arb(; prec), Arb(x; prec); prec)
    lΓy = Arblib.lgamma!(Arb(; prec), Arb(y; prec); prec)
    lΓxy = Arblib.lgamma!(Arb(; prec), Arb(x + y; prec); prec)
    return T(exp(lΓx + lΓy - lΓxy))
end

@doc raw"""
Volume of a super-spheroid:

```math
V = 4\,a^2 c\, n^2 \, B\!\left(\tfrac{n}{2}+1,\, n\right) B\!\left(\tfrac{n}{2},\, \tfrac{n+2}{2}\right)
```

where ``B`` is the beta function.  At ``n=1`` this reduces to ``\tfrac{4}{3}\pi a^2 c``
(spheroid).
"""
function volume(s::SuperSpheroid{T}) where {T}
    n = s.n
    b1 = _beta_lgamma(n / 2 + 1, n)
    b2 = _beta_lgamma(n / 2, (n + 2) / 2)
    return 4 * s.a^2 * s.c * n^2 * b1 * b2
end

volume_equivalent_radius(s::SuperSpheroid) = ∛(3 * volume(s) / (4 * oftype(s.a, π)))

has_symmetric_plane(::SuperSpheroid) = true

@testitem "SuperSpheroid utility functions" begin
    using TransitionMatrices: SuperSpheroid, volume, volume_equivalent_radius,
                              has_symmetric_plane, Double64, Float128, ComplexF128

    # n=1 reduces to spheroid: V = 4π/3 * a^2 * c
    @testset "n=1 reduction to spheroid" begin
        a, c = 1.0, 2.0
        s = SuperSpheroid(a, c, 1.0, complex(1.5))
        @test volume(s) ≈ 4π / 3 * a^2 * c rtol=1e-8
        @test volume_equivalent_radius(s) ≈ ∛(a^2 * c) rtol=1e-8
    end

    # `volume` keeps the shape's real type and its full precision for every supported
    # numeric type (regression guard for the Arb→T beta computation: no BigFloat detour,
    # instance precision honoured, and the Quadmath `precision(::Float128)` workaround).
    @testset "volume type & precision for $T" for (T, CT, mk) in (
        (Float64, ComplexF64, 1e-14),
        (Double64, Complex{Double64}, 1e-28),
        (Float128, ComplexF128, 1e-30),
        (BigFloat, Complex{BigFloat}, 1e-70))
        a, c = T(12) / 10, T(8) / 10
        s = SuperSpheroid(a, c, one(T), complex(CT(3) / 2))
        v = volume(s)
        @test v isa T
        @test abs(v - 4 * T(π) / 3 * a^2 * c) < mk * abs(v)
    end

    # n=2 reduces to octahedron: V = 4/3 * a^2 * c
    @testset "n=2 reduction to octahedron" begin
        a, c = 1.0, 1.0
        s = SuperSpheroid(a, c, 2.0, complex(1.5))
        @test volume(s) ≈ 4 / 3 * a^2 * c rtol=1e-8
    end

    # Volume vs Monte-Carlo ∈ test
    @testset "Volume vs Monte-Carlo for n=$n" for n in [1.2, 1.7, 2.5]
        a, c = 1.2, 0.8
        s = SuperSpheroid(a, c, n, complex(1.5))
        V_analytic = volume(s)

        rng_state = 42
        N = 200_000
        # Sample within bounding box [-a,a]^2 x [-c,c]
        hits = 0
        let seed = UInt64(rng_state)
            # Simple LCG for reproducibility
            x = zeros(3)
            for _ in 1:N
                seed = seed * 6364136223846793005 + 1442695040888963407
                x[1] = (Float64(seed >> 11) / 2^53 * 2 - 1) * a
                seed = seed * 6364136223846793005 + 1442695040888963407
                x[2] = (Float64(seed >> 11) / 2^53 * 2 - 1) * a
                seed = seed * 6364136223846793005 + 1442695040888963407
                x[3] = (Float64(seed >> 11) / 2^53 * 2 - 1) * c
                hits += x ∈ s ? 1 : 0
            end
        end
        V_box = 2a * 2a * 2c
        V_mc = V_box * hits / N
        @test abs(V_analytic - V_mc) / V_analytic < 0.01
    end

    @test has_symmetric_plane(SuperSpheroid(1.0, 1.0, 1.5, complex(1.5)))
end

function rmax(s::SuperSpheroid)
    # Equatorial max is a (at cardinal phi=0,pi/2), polar max is c.
    # For n>=1 the equatorial maximum is a (super-circle corners are at axes).
    # For n<1 the diagonal is larger, but IITM convergence typically uses n>=1.
    n = Float64(s.n)
    if n >= 1
        return max(s.a, s.c)
    else
        # equatorial diagonal: a * 2^((1-n)/2)
        return max(s.a * oftype(s.a, 2)^((1 - n) / 2), s.c)
    end
end

function rmin(s::SuperSpheroid)
    # Inscribed-sphere radius = min over directions û of the surface radius r(û). There is
    # NO simple closed form: for concave superquadrics (n>2) with a≠c the closest surface
    # point sits at a transcendental, off-symmetry direction (the symmetry-direction minimum
    # over-estimates by a few %). So minimise r(û) over a dense direction grid. An over-large
    # rₘᵢₙ — e.g. min(a,c), which ignores the body-diagonal/face dimples — makes the IITM
    # Mie-init sphere poke outside the particle and gives wrong cross sections, so a small
    # shrink (below) keeps it safely inscribed.
    n = Float64(s.n)
    a = Float64(s.a)
    c = Float64(s.c)
    q = 2 / n
    M = 90
    rmn = Inf
    for i in 0:M, j in 0:M
        ϑ = i * (π / 2) / M
        φ = j * (π / 2) / M
        ux = sin(ϑ) * cos(φ)
        uy = sin(ϑ) * sin(φ)
        uz = cos(ϑ)
        g = (ux / a)^q + (uy / a)^q + (uz / c)^q
        rmn = min(rmn, g^(-n / 2))
    end
    # The grid minimum can only OVER-estimate the true minimum (it may miss the exact
    # closest direction between nodes, by ≲1e-3 here). rₘᵢₙ must never exceed the true
    # inscribed radius — an over-large value pokes the Mie-init sphere outside the
    # particle — so shrink by a margin comfortably larger than the grid error.
    return oftype(s.a, 0.99 * rmn)
end

function Base.:∈(x, s::SuperSpheroid)
    return (abs(x[1]) / s.a)^(2 / s.n) +
           (abs(x[2]) / s.a)^(2 / s.n) +
           (abs(x[3]) / s.c)^(2 / s.n) <= 1
end

refractive_index(s::SuperSpheroid, x) = x ∈ s ? s.m : one(s.m)

@testitem "SuperSpheroid n=1 reduces to a spheroid (IITM)" begin
    using TransitionMatrices: SuperSpheroid, Spheroid, transition_matrix, IITM,
                              calc_Csca, calc_Cext

    m = complex(1.5)
    # n=1 ⇒ the super-spheroid IS a spheroid. Compare both through the SAME IITM solver
    # with a matched explicit rₘᵢₙ: this checks the geometric reduction directly (to
    # machine precision) and is robust to the rₘᵢₙ default and to cross-method
    # (EBCM-vs-IITM) discretisation differences.
    slv = IITM(5, 60, 200, 100; rₘᵢₙ = 0.6)
    𝐓sph = transition_matrix(Spheroid(1.0, 2.0, m), 2π, slv)
    𝐓ss = transition_matrix(SuperSpheroid(1.0, 2.0, 1.0, m), 2π, slv)

    @test abs(calc_Csca(𝐓sph) - calc_Csca(𝐓ss)) < 1e-10
    @test abs(calc_Cext(𝐓sph) - calc_Cext(𝐓ss)) < 1e-10
end

@testitem "SuperSpheroid geometry predicates (incl. n<1 rmax branch)" begin
    using TransitionMatrices: SuperSpheroid
    for n in (0.7, 1.0, 2.5)   # n<1 hits the diagonal-extended rmax branch
        s = SuperSpheroid(1.0, 2.0, n, complex(1.5))
        @test TransitionMatrices.rmax(s) ≥ TransitionMatrices.rmin(s) > 0
        @test [0.0, 0.0, 0.0] ∈ s
        @test [5.0, 0.0, 0.0] ∉ s
        @test TransitionMatrices.refractive_index(s, [0.0, 0.0, 0.0]) == complex(1.5)
        @test TransitionMatrices.refractive_index(s, [5.0, 0.0, 0.0]) == one(complex(1.5))
    end
end
