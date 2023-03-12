@doc raw"""
A spheroidal scatterer.

Attributes:

- `a`: the length of the semi-major axis.
- `c`: the length of the semi-minor axis.
- `m`: the relative complex refractive index.
"""
struct Spheroid{T, CT} <: AbstractAxisymmetricShape{T, CT}
    a::T
    c::T
    m::CT
end

volume(s::Spheroid) = 4.0 / 3.0 * π * s.a^2 * s.c
volume_equivalent_radius(s::Spheroid) = ∛(s.a^2 * s.c)
has_symmetric_plane(::Spheroid) = true

@testitem "Utility functions are correct" begin
    using TransitionMatrices: Spheroid, volume, volume_equivalent_radius,
                              has_symmetric_plane

    @testset "a = $a, c = $c" for (a, c, V, rᵥ) in [
        (1.0, 0.5, 2.0943951023931953, 0.7937005259840998),
        (1.0, 1.0, 4.1887902047863905, 1.0),
        (0.5, 1.0, 1.0471975511965976, 0.6299605249474366),
    ]
        s = Spheroid(a, c, 1.311)
        @test volume(s) ≈ V
        @test volume_equivalent_radius(s) ≈ rᵥ
        @test has_symmetric_plane(s)
    end
end

@doc raw"""
```
gaussquad(s::Spheroid{T}, ngauss) where {T}
```

Evaluate the quadrature points, weights and the corresponding radius and radius derivative (to ``\vartheta``) for a spheroid.

Returns: (`x`, `w`, `r`, `r′`)

- `x`: the quadrature points.
- `w`: the quadrature weights.
- `r`: the radius at each quadrature point.
- `r′`: the radius derivative at each quadrature point.
"""
function gaussquad(s::Spheroid{T}, ngauss) where {T}
    x, w = gausslegendre(T, ngauss)
    r = similar(x)
    r′ = similar(x)

    # We could have used @turbo for primitive types like 
    # `Float64`, but this is not the bottleneck, so we only use 
    # `@simd` here.
    @simd for i in 1:(ngauss ÷ 2)
        cosϑ = x[i]
        sinϑ = √(1.0 - cosϑ^2)
        r[i] = s.a * s.c * √(1.0 / (s.a^2 * cosϑ^2 + s.c^2 * sinϑ^2))
        r[ngauss + 1 - i] = r[i]
        r′[i] = r[i]^3 * cosϑ * sinϑ * (s.a^2 - s.c^2) / (s.a^2 * s.c^2)
        r′[ngauss + 1 - i] = -r′[i]
    end

    return x, w, r, r′
end

@testitem "radius derivative has correct sign" begin
    using TransitionMatrices: Spheroid, gaussquad
    using Arblib: Arb

    @testset "a = 1.0, c = 2.0, T = $T" for T in (Float64, Arb)
        s = Spheroid(T(1.0), T(2.0), T(1.0))
        _, _, _, r′ = gaussquad(s, 100)

        # r′[1] is near pi (cosϑ ≈ -1)
        @test r′[1] > 0.0

        # r′[end] is near 0 (cosϑ ≈ 1)
        @test r′[end] < 0.0
    end

    @testset "a = 2.0, c = 1.0, T = $T" for T in (Float64, Arb)
        s = Spheroid(T(2.0), T(1.0), T(1.0))
        _, _, _, r′ = gaussquad(s, 100)

        # r′[1] is near pi (cosϑ ≈ -1)
        @test r′[1] < 0.0

        # r′[end] is near 0 (cosϑ ≈ 1)
        @test r′[end] > 0.0
    end

    @testset "a = 1.0, c = 1.0, T = $T" for T in (Float64, Arb)
        s = Spheroid(T(1.0), T(1.0), T(1.0))
        _, _, r, r′ = gaussquad(s, 100)

        # Since this is actually a sphere, r should be equal everywhere, while r′ should be zero.
        @test all(r .≈ 1.0)
        @test all(iszero.(r′))
    end
end

function estimate_integral_loss(s::Spheroid, nₘₐₓ)
    mach = machine("fixtures/spheroids.Q.model")
    predict(mach,
            [(rmin = min(s.a, s.c), rmax = max(s.a, s.c), mr = real(s.m), mi = imag(s.m),
              nmax = nₘₐₓ)])[1]
end

function estimate_inverse_loss(s::Spheroid, nₘₐₓ)
    mach = machine("fixtures/spheroids.T.model")
    predict(mach,
            [(rmin = min(s.a, s.c), rmax = max(s.a, s.c), mr = real(s.m), mi = imag(s.m),
              nmax = nₘₐₓ)])[1]
end

function estimate_total_loss(s::Spheroid, nₘₐₓ)
    mach = machine("fixtures/spheroids.AT.model")
    predict(mach,
            [(rmin = min(s.a, s.c), rmax = max(s.a, s.c), mr = real(s.m), mi = imag(s.m),
              nmax = nₘₐₓ)])[1]
end
