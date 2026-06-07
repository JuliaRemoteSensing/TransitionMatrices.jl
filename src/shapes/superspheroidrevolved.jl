@doc raw"""
An axisymmetric super-spheroid scatterer (surface of revolution of a superellipse).

This is the surface of revolution

```math
\left(\frac{\varrho}{a}\right)^{2/n} + \left(\frac{z}{c}\right)^{2/n} = 1,
\qquad \varrho = \sqrt{x^2 + y^2}
```

Unlike `SuperSpheroid` (which has 4-fold D4h symmetry from super-circular cross-sections),
this shape is **truly axisymmetric** and is compatible with the EBCM / Sh-matrix solvers.

At ``n = 1`` it reduces exactly to `Spheroid(a, c)`.

Conventions:
- ``n = 1``: ordinary spheroid.
- ``n = 2``: bicone (diamond/spindle shape).
- ``1 < n < 2``: convex rounded.
- ``n > 2``: concave (waisted), analogous to the Lin/Bi dust regime.

Attributes:

- `a`: the equatorial semi-axis.
- `c`: the polar semi-axis.
- `n`: the roundness exponent.
- `m`: the relative complex refractive index.
"""
struct SuperSpheroidRevolved{T, CT} <: AbstractAxisymmetricShape{T, CT}
    a::T
    c::T
    n::T
    m::CT
end

SuperSpheroidRevolved(a, c, n, m) = SuperSpheroidRevolved{typeof(a), typeof(m)}(a, c, n, m)

@doc raw"""
Volume of the axisymmetric super-spheroid (surface of revolution):

```math
V = \pi a^2 c \, n \, B\!\left(\tfrac{n}{2},\, n+1\right)
```

At ``n = 1``: ``B(1/2, 2) = 4/3``, so ``V = \tfrac{4}{3}\pi a^2 c`` (spheroid).
At ``n = 2``: ``B(1, 3) = 1/3``, so ``V = \tfrac{2}{3}\pi a^2 c`` (bicone).
"""
function volume(s::SuperSpheroidRevolved{T}) where {T}
    n = s.n
    b = _beta_lgamma(n / 2, n + 1)
    return oftype(s.a, π) * s.a^2 * s.c * n * b
end

volume_equivalent_radius(s::SuperSpheroidRevolved) = ∛(3 * volume(s) / (4 * oftype(s.a, π)))

has_symmetric_plane(::SuperSpheroidRevolved) = true

@testitem "SuperSpheroidRevolved utility functions" begin
    using TransitionMatrices: SuperSpheroidRevolved, Spheroid, volume,
                              volume_equivalent_radius, has_symmetric_plane

    # n=1 must exactly match Spheroid volume
    @testset "n=1 reduction to spheroid" begin
        a, c = 1.0, 2.0
        s = SuperSpheroidRevolved(a, c, 1.0, complex(1.5))
        sph = Spheroid(a, c, complex(1.5))
        @test volume(s) ≈ volume(sph) rtol=1e-8
        @test volume_equivalent_radius(s) ≈ volume_equivalent_radius(sph) rtol=1e-8
    end

    # n=2 → bicone: V = 2π/3 * a^2 * c
    @testset "n=2 reduction to bicone" begin
        a, c = 1.0, 1.5
        s = SuperSpheroidRevolved(a, c, 2.0, complex(1.5))
        @test volume(s) ≈ 2π / 3 * a^2 * c rtol=1e-8
    end

    # Volume vs numeric integration (Simpson on the revolved surface)
    @testset "Volume vs quadrature for n=$n" for n in [1.2, 1.7, 2.5]
        a, c = 1.2, 0.8
        s = SuperSpheroidRevolved(a, c, n, complex(1.5))
        V_analytic = volume(s)

        # V = π ∫_{-c}^{c} ϱ(z)^2 dz  where (ϱ/a)^(2/n) + (z/c)^(2/n) = 1
        # ϱ(z) = a * (1 - (|z|/c)^(2/n))^(n/2)
        Npts = 10_000
        dz = 2c / Npts
        V_num = 0.0
        for i in 0:(Npts - 1)
            z = -c + (i + 0.5) * dz
            t = abs(z) / c
            t < 1 || continue
            ϱ = a * (1 - t^(2 / n))^(n / 2)
            V_num += π * ϱ^2 * dz
        end
        @test abs(V_analytic - V_num) / V_analytic < 1e-3
    end

    @test has_symmetric_plane(SuperSpheroidRevolved(1.0, 1.0, 1.5, complex(1.5)))
end

@doc raw"""
```
gaussquad(s::SuperSpheroidRevolved{T}, ngauss) where {T}
```

Evaluate the quadrature points, weights and the corresponding radius ``r(\vartheta)``
and its derivative ``r'(\vartheta) = \mathrm{d}r/\mathrm{d}\vartheta`` for an
axisymmetric super-spheroid.

The surface in spherical coordinates is defined by ``r(\vartheta)`` satisfying

```math
\left(\frac{r\sin\vartheta}{a}\right)^{2/n} + \left(\frac{r\cos\vartheta}{c}\right)^{2/n} = 1
```

Let ``p = 2/n`` and ``g(\vartheta) = (\sin\vartheta/a)^p + (\cos\vartheta/c)^p``.  Then

```math
r(\vartheta) = g(\vartheta)^{-1/p}, \qquad
r'(\vartheta) = -g^{-1/p-1}\left[\frac{\sin^{p-1}\vartheta\cos\vartheta}{a^p}
                                   - \frac{\cos^{p-1}\vartheta\sin\vartheta}{c^p}\right]
```

Returns ``(x, w, r, r')`` where ``x = \cos\vartheta`` are Gauss–Legendre nodes on ``[-1,1]``.
"""
function gaussquad(s::SuperSpheroidRevolved{T}, ngauss) where {T}
    x, w = gausslegendre(T, ngauss)
    r = similar(x)
    r′ = similar(x)

    p = 2 / s.n   # exponent
    ap = s.a^p
    cp = s.c^p

    @simd for i in 1:ngauss
        cosϑ = x[i]
        sinϑ = √(1 - cosϑ^2)

        # g(ϑ) = (sinϑ/a)^p + (cosϑ/c)^p  (using |cosϑ|: pole symmetry handled below)
        # Use absolute values for generality; z-symmetry from has_symmetric_plane
        sinp = abs(sinϑ)^p / ap
        cosp = abs(cosϑ)^p / cp
        g = sinp + cosp

        # r = g^(-1/p) = g^(-n/2)
        r[i] = g^(-s.n / 2)

        # r′ = -g^(-1/p - 1) * [sinϑ^(p-1)*cosϑ/a^p - cosϑ^(p-1)*sinϑ/c^p]
        # = -g^(-n/2 - 1) * [sign(cosϑ)*|sinϑ|^(p-1)*|cosϑ|^1/a^p
        #                     - sign(sinϑ)*|cosϑ|^(p-1)*|sinϑ|^1/c^p]
        # For ϑ ∈ (0,π): sinϑ ≥ 0, sign(sinϑ)=+1
        # cosϑ can be negative (ϑ > π/2)
        if sinϑ == 0
            # Defensive pole guard: gausslegendre nodes are interior (|cosϑ|<1), so
            # sinϑ>0 always and this branch is unreachable via `gaussquad` — excluded
            # from coverage rather than removed (keeps r′ finite if ever called at a pole).
            r′[i] = zero(T)  # COV_EXCL_LINE
        else
            # sinϑ > 0 for i not at exactly the pole
            dg_dϑ = (p * abs(sinϑ)^(p - 1) * cosϑ / ap -
                     p * abs(cosϑ)^(p - 1) * sinϑ * sign(cosϑ) / cp)
            # r′ = dr/dϑ = -1/p * g^(-1/p-1) * dg/dϑ = -(n/2) * g^(-n/2-1) * dg/dϑ
            r′[i] = -(s.n / 2) * g^(-s.n / 2 - 1) * dg_dϑ
        end
    end

    return x, w, r, r′
end

@testitem "SuperSpheroidRevolved gaussquad matches Spheroid at n=1" begin
    using TransitionMatrices: SuperSpheroidRevolved, Spheroid, gaussquad

    a, c = 1.0, 2.0
    s = SuperSpheroidRevolved(a, c, 1.0, complex(1.5))
    sph = Spheroid(a, c, complex(1.5))
    ngauss = 60

    _, _, r_s, rp_s = gaussquad(s, ngauss)
    _, _, r_sph, rp_sph = gaussquad(sph, ngauss)

    @test all(r_s .≈ r_sph)
    @test all(rp_s .≈ rp_sph)
end

@testitem "SuperSpheroidRevolved gaussquad r' sign" begin
    using TransitionMatrices: SuperSpheroidRevolved, gaussquad

    # Prolate (c > a): at theta=pi (costheta=-1, south pole), dr/dtheta > 0
    # (r decreases from c to a as theta goes from pi to pi/2)
    @testset "prolate a=$a c=$c n=$n" for (a, c, n) in [(1.0, 2.0, 1.0), (1.0, 2.0, 1.5)]
        s = SuperSpheroidRevolved(a, c, n, complex(1.5))
        _, _, _, r′ = gaussquad(s, 100)
        @test r′[1] > 0    # near theta=pi (costheta ≈ -1)
        @test r′[end] < 0  # near theta=0  (costheta ≈ +1)
    end
end

function rmax(s::SuperSpheroidRevolved)
    # r(ϑ) = g(ϑ)^(-n/2); extremes occur at poles (ϑ=0,π) and equator (ϑ=π/2)
    # r at pole: g = (1/c)^p → r = c
    # r at equator: g = (1/a)^p → r = a
    # For n>2 (concave), rmax may occur at an interior ϑ; compute numerically.
    return _superspheroidrevolved_rmax_rmin(s)[1]
end

function rmin(s::SuperSpheroidRevolved)
    return _superspheroidrevolved_rmax_rmin(s)[2]
end

function _superspheroidrevolved_rmax_rmin(s::SuperSpheroidRevolved)
    p = 2 / s.n
    ap = s.a^p
    cp = s.c^p
    # Sample r over ϑ ∈ [0, π/2] (symmetric about equator and pole)
    N = 1000
    rvals = [begin
        ϑ = i * π / 2 / N
        sinϑ = sin(ϑ)
        cosϑ = cos(ϑ)
        g = sinϑ^p / ap + cosϑ^p / cp
        g^(-s.n / 2)
    end
             for i in 0:N]
    return maximum(rvals), minimum(rvals)
end

function Base.:∈(x, s::SuperSpheroidRevolved)
    ϱ = √(x[1]^2 + x[2]^2)
    return (ϱ / s.a)^(2 / s.n) + (abs(x[3]) / s.c)^(2 / s.n) <= 1
end

refractive_index(s::SuperSpheroidRevolved, x) = x ∈ s ? s.m : one(s.m)

@testitem "SuperSpheroidRevolved n=1 EBCM T-matrix matches Spheroid" begin
    using TransitionMatrices: SuperSpheroidRevolved, Spheroid, EBCM,
                              transition_matrix, calc_Csca, calc_Cext

    a, c = 1.0, 2.0
    m = complex(1.5)

    sph = Spheroid(a, c, m)
    ssr = SuperSpheroidRevolved(a, c, 1.0, m)

    nₘₐₓ = 5
    Ng = 80

    𝐓sph = transition_matrix(sph, 2π, EBCM(nₘₐₓ, Ng))
    𝐓ssr = transition_matrix(ssr, 2π, EBCM(nₘₐₓ, Ng))

    @test abs(calc_Csca(𝐓sph) - calc_Csca(𝐓ssr)) < 1e-10
    @test abs(calc_Cext(𝐓sph) - calc_Cext(𝐓ssr)) < 1e-10
end

@testitem "SuperSpheroidRevolved geometry predicates (∈, rmax, rmin, refractive_index)" begin
    using TransitionMatrices: SuperSpheroidRevolved
    # span n<1, ellipse, convex, concave so rmax/rmin sampling + ∈ are exercised
    for n in (0.8, 1.0, 1.5, 2.5)
        s = SuperSpheroidRevolved(1.0, 2.0, n, complex(1.5))
        rmx = TransitionMatrices.rmax(s)
        rmn = TransitionMatrices.rmin(s)
        @test rmx ≥ rmn > 0
        @test [0.0, 0.0, 0.0] ∈ s          # centre inside
        @test [5.0, 0.0, 0.0] ∉ s          # far outside
        @test [0.0, 0.0, 0.99 * s.c] ∈ s   # near a pole, inside
        @test TransitionMatrices.refractive_index(s, [0.0, 0.0, 0.0]) == complex(1.5)
        @test TransitionMatrices.refractive_index(s, [5.0, 0.0, 0.0]) == one(complex(1.5))
    end
end
