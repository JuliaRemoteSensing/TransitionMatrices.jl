@doc raw"""
A super-ellipsoid scatterer (Bi et al. 2018 D2h symmetry).

Defined by the implicit surface

```math
\left[\left|\frac{x}{a}\right|^{2/e} + \left|\frac{y}{b}\right|^{2/e}\right]^{e/n}
+ \left|\frac{z}{c}\right|^{2/n} = 1
```

For ``a = b`` this reduces to a 4-fold (D4h) shape; for ``a \neq b`` it is D2h (2-fold).
At ``e = n = 1`` it reduces to an ordinary ellipsoid.

Conventions (Bi et al. 2018 Opt. Express, doi:10.1364/OE.26.001726):
- ``e``: east–west (equatorial cross-section) roundness exponent.
- ``n``: north–south (meridional profile) roundness exponent.
- ``e = n = 1``: ordinary ellipsoid.
- ``e = n = 2``: octahedron-like shape.

Attributes:

- `a`: the ``x`` semi-axis.
- `b`: the ``y`` semi-axis.
- `c`: the ``z`` (polar) semi-axis.
- `e`: the equatorial roundness exponent.
- `n`: the meridional roundness exponent.
- `m`: the relative complex refractive index.
"""
struct SuperEllipsoid{T, CT} <: AbstractNFoldShape{2, T, CT}
    a::T
    b::T
    c::T
    e::T
    n::T
    m::CT
end

SuperEllipsoid(a, b, c, e, n, m) = SuperEllipsoid{typeof(a), typeof(m)}(a, b, c, e, n, m)

# Helper shared across super-spheroid shapes (defined in superspheroid.jl):
# _beta_lgamma(x, y) — uses Arblib lgamma!

@doc raw"""
Volume of a super-ellipsoid:

```math
V = a\,b\,c\, e\,n \, B\!\left(\tfrac{e}{2},\, \tfrac{e}{2}\right)
                       B\!\left(\tfrac{n}{2},\, n+1\right)
```

where ``B(x,y) = \Gamma(x)\Gamma(y)/\Gamma(x+y)`` is the beta function.
At ``e = n = 1`` this reduces to ``\tfrac{4}{3}\pi abc`` (ellipsoid).

Derivation: the cross-section at height ``z`` is the super-ellipse
``\{|x/a|^{2/e}+|y/b|^{2/e} \le (1-|z/c|^{2/n})^{n/e}\}``
whose area is ``ab\,(1-|z/c|^{2/n})^n \cdot e\,B(e/2,e/2)``;
integrating gives the formula above.
"""
function volume(s::SuperEllipsoid{T}) where {T}
    e = s.e
    n = s.n
    b1 = _beta_lgamma(e / 2, e / 2)   # B(e/2, e/2) = Γ(e/2)²/Γ(e)
    b2 = _beta_lgamma(n / 2, n + 1)   # B(n/2, n+1)
    return s.a * s.b * s.c * e * n * b1 * b2
end

volume_equivalent_radius(s::SuperEllipsoid) = ∛(3 * volume(s) / (4 * oftype(s.a, π)))

has_symmetric_plane(::SuperEllipsoid) = true

@testitem "SuperEllipsoid utility functions" begin
    using TransitionMatrices: SuperEllipsoid, volume, volume_equivalent_radius,
                              has_symmetric_plane

    # e=n=1 reduces to ellipsoid: V = 4π/3 * a * b * c
    @testset "e=n=1 reduction to ellipsoid" begin
        a, b, c = 1.0, 1.5, 2.0
        s = SuperEllipsoid(a, b, c, 1.0, 1.0, complex(1.5))
        @test volume(s) ≈ 4π / 3 * a * b * c rtol=1e-8
        @test volume_equivalent_radius(s) ≈ ∛(a * b * c) rtol=1e-8
    end

    # a=b, e=n=1 reduces to spheroid
    @testset "a=b, e=n=1 reduces to spheroid" begin
        a, c = 1.0, 2.0
        s = SuperEllipsoid(a, a, c, 1.0, 1.0, complex(1.5))
        @test volume(s) ≈ 4π / 3 * a^2 * c rtol=1e-8
    end

    # Volume vs independent numerical integration via cross-sectional slices
    # V = ∫_{-c}^{c} A(z) dz, A(z) = ab·R(z)^e·A_unit
    # R(z) = (1-(|z|/c)^(2/n))^(n/e), A_unit = area({|u|^(2/e)+|v|^(2/e)≤1}) = 4·∫_0^1(1-u^(2/e))^(e/2)du
    @testset "Volume vs cross-section integration for e=$e, n=$n" for (e, n) in [(1.5, 1.2), (2.0, 2.0)]
        a, b, c = 1.0, 1.2, 0.8
        s = SuperEllipsoid(a, b, c, e, n, complex(1.5))
        V_analytic = volume(s)

        Npts = 5000
        du = 1.0 / Npts
        A_unit = 4 * sum((1 - ((i - 0.5) * du)^(2 / e))^(e / 2) for i in 1:Npts) * du

        dz = 2c / Npts
        V_num = 0.0
        for i in 1:Npts
            z = -c + (i - 0.5) * dz
            t = abs(z) / c
            t < 1 || continue
            R = (1 - t^(2 / n))^(n / e)
            V_num += a * b * R^e * A_unit * dz
        end

        @test abs(V_analytic - V_num) / V_analytic < 0.005
    end

    # Fully independent 3D Monte-Carlo via `∈` (shares no derivation with the closed-form
    # beta-function volume — guards the corrected SuperEllipsoid formula). Deterministic
    # LCG seed so the test is reproducible, not flaky.
    @testset "Volume vs 3D Monte-Carlo for e=$e, n=$n" for (e, n) in [(0.8, 1.6), (1.5, 1.2), (2.0, 2.5)]
        a, b, c = 1.0, 1.2, 0.8
        s = SuperEllipsoid(a, b, c, e, n, complex(1.5))
        V_analytic = volume(s)
        N = 400_000
        hits = 0
        let seed = UInt64(20240607)
            x = zeros(3)
            for _ in 1:N
                seed = seed * 6364136223846793005 + 1442695040888963407
                x[1] = (Float64(seed >> 11) / 2^53 * 2 - 1) * a
                seed = seed * 6364136223846793005 + 1442695040888963407
                x[2] = (Float64(seed >> 11) / 2^53 * 2 - 1) * b
                seed = seed * 6364136223846793005 + 1442695040888963407
                x[3] = (Float64(seed >> 11) / 2^53 * 2 - 1) * c
                hits += x ∈ s ? 1 : 0
            end
        end
        V_mc = 2a * 2b * 2c * hits / N
        @test abs(V_analytic - V_mc) / V_analytic < 0.01
    end

    @test has_symmetric_plane(SuperEllipsoid(1.0, 1.2, 0.8, 1.5, 1.2, complex(1.5)))
end

function rmax(s::SuperEllipsoid)
    # For e,n ≥ 1: the surface never extends beyond the cardinal semi-axes,
    # so rmax = max(a, b, c).
    # For e or n < 1 (cube-like bumpy shapes): corners can extend beyond axes;
    # hypot gives a safe conservative upper bound in that case.
    if s.e >= 1 && s.n >= 1
        return max(s.a, s.b, s.c)
    else
        return hypot(s.a, s.b, s.c)
    end
end

function rmin(s::SuperEllipsoid)
    # Inscribed-sphere radius = min over directions û of the surface radius
    #   r(û) = ((|ux/a|^(2/e)+|uy/b|^(2/e))^(e/n) + |uz/c|^(2/n))^(-n/2).
    # There is NO simple closed form: for concave superquadrics the closest surface point
    # sits at a transcendental, off-symmetry direction (the symmetry-direction minimum
    # over-estimates by a few % for e≠n). So minimise r(û) over a dense direction grid.
    # An over-large rₘᵢₙ — e.g. min(a,b,c), which ignores the body-diagonal/face dimples —
    # makes the IITM Mie-init sphere poke outside the particle and gives wrong cross
    # sections, so a small shrink (below) keeps it safely inscribed.
    e = Float64(s.e)
    n = Float64(s.n)
    a = Float64(s.a)
    b = Float64(s.b)
    c = Float64(s.c)
    M = 90
    rmn = Inf
    for i in 0:M, j in 0:M
        ϑ = i * (π / 2) / M
        φ = j * (π / 2) / M
        ux = sin(ϑ) * cos(φ)
        uy = sin(ϑ) * sin(φ)
        uz = cos(ϑ)
        Σ = (ux / a)^(2 / e) + (uy / b)^(2 / e)
        g = Σ^(e / n) + (uz / c)^(2 / n)
        rmn = min(rmn, g^(-n / 2))
    end
    # The grid minimum can only OVER-estimate the true minimum (it may miss the exact
    # closest direction between nodes). rₘᵢₙ must never exceed the true inscribed radius —
    # an over-large value pokes the Mie-init sphere outside the particle — so shrink by a
    # margin comfortably larger than the grid error.
    return oftype(s.a, 0.99 * rmn)
end

function Base.:∈(x, s::SuperEllipsoid)
    inner = (abs(x[1]) / s.a)^(2 / s.e) + (abs(x[2]) / s.b)^(2 / s.e)
    return inner^(s.e / s.n) + (abs(x[3]) / s.c)^(2 / s.n) <= 1
end

refractive_index(s::SuperEllipsoid, x) = x ∈ s ? s.m : one(s.m)

@testitem "SuperEllipsoid e=n=1 (a=b) reduces to a spheroid (IITM)" begin
    using TransitionMatrices: SuperEllipsoid, Spheroid, transition_matrix, IITM,
                              calc_Csca, calc_Cext

    m = complex(1.5)
    # a=b, e=n=1 ⇒ the super-ellipsoid IS a spheroid. Compare both through the SAME IITM
    # solver with a matched explicit rₘᵢₙ: tests the geometric reduction directly (to
    # machine precision), robust to the rₘᵢₙ default and cross-method discretisation.
    slv = IITM(5, 60, 200, 100; rₘᵢₙ = 0.6)
    𝐓sph = transition_matrix(Spheroid(1.0, 2.0, m), 2π, slv)
    𝐓se = transition_matrix(SuperEllipsoid(1.0, 1.0, 2.0, 1.0, 1.0, m), 2π, slv)

    @test abs(calc_Csca(𝐓sph) - calc_Csca(𝐓se)) < 1e-10
    @test abs(calc_Cext(𝐓sph) - calc_Cext(𝐓se)) < 1e-10
end

@testitem "SuperEllipsoid geometry predicates (incl. e,n<1 rmax branch)" begin
    using TransitionMatrices: SuperEllipsoid
    for (e, n) in ((0.7, 1.5), (1.0, 1.0), (2.0, 2.5))   # e<1 hits the hypot rmax branch
        s = SuperEllipsoid(1.0, 1.2, 0.8, e, n, complex(1.5))
        @test TransitionMatrices.rmax(s) ≥ TransitionMatrices.rmin(s) > 0
        @test [0.0, 0.0, 0.0] ∈ s
        @test [5.0, 0.0, 0.0] ∉ s
        @test TransitionMatrices.refractive_index(s, [0.0, 0.0, 0.0]) == complex(1.5)
        @test TransitionMatrices.refractive_index(s, [5.0, 0.0, 0.0]) == one(complex(1.5))
    end
end
