export Spheroid

@doc raw"""
A spheroid scatterer.

Attributes:

- `a`: Length of the semi-major axis
- `c`: Length of the semi-minor axis
- `m`: Relative complex refractive index
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

function radius_and_deriv!(r, dr, s::Spheroid{Float64}, x)
    ngauss = length(x)

    # Use @turbo for Float64 only (loop vectorization does not work well for other types, e.g., Double64, Float128, and does not work at all for Arb)
    @turbo for i in 1:(ngauss ÷ 2)
        cosϑ = x[i]
        sinϑ = √(1.0 - cosϑ^2)
        r[i] = s.a * s.c * √(1.0 / (s.a^2 * cosϑ^2 + s.c^2 * sinϑ^2))
        r[ngauss + 1 - i] = r[i]
        dr[i] = r[i]^3 * cosϑ * sinϑ * (s.a^2 - s.c^2) / (s.a^2 * s.c^2)
        dr[ngauss + 1 - i] = -dr[i]
    end
end

function radius_and_deriv!(r, dr, s::Spheroid{T}, x) where {T}
    ngauss = length(x)

    # For other types, we only use @inbounds to skip bounds checking
    @inbounds for i in 1:(ngauss ÷ 2)
        cosϑ = x[i]
        sinϑ = √(1.0 - cosϑ^2)
        r[i] = s.a * s.c * √(1.0 / (s.a^2 * cosϑ^2 + s.c^2 * sinϑ^2))
        r[ngauss + 1 - i] = r[i]
        dr[i] = r[i]^3 * cosϑ * sinϑ * (s.a^2 - s.c^2) / (s.a^2 * s.c^2)
        dr[ngauss + 1 - i] = -dr[i]
    end
end

@testitem "radius derivative has correct sign" begin
    using TransitionMatrices: Spheroid, gausslegendre, radius_and_deriv!
    using Arblib: Arb

    @testset "a = 1.0, c = 2.0, T = $T" for T in (Float64, Arb)
        s = Spheroid(T(1.0), T(2.0), T(1.0))
        x, w = gausslegendre(T, 100)
        r = similar(x)
        dr = similar(x)
        radius_and_deriv!(r, dr, s, x)

        # dr[1] is near pi (cosϑ ≈ -1)
        @test dr[1] > 0.0

        # dr[end] is near 0 (cosϑ ≈ 1)
        @test dr[end] < 0.0
    end

    @testset "a = 2.0, c = 1.0, T = $T" for T in (Float64, Arb)
        s = Spheroid(T(2.0), T(1.0), T(1.0))
        x, w = gausslegendre(T, 100)
        r = similar(x)
        dr = similar(x)
        radius_and_deriv!(r, dr, s, x)

        # dr[1] is near pi (cosϑ ≈ -1)
        @test dr[1] < 0.0

        # dr[end] is near 0 (cosϑ ≈ 1)
        @test dr[end] > 0.0
    end

    @testset "a = 1.0, c = 1.0, T = $T" for T in (Float64, Arb)
        s = Spheroid(T(1.0), T(1.0), T(1.0))
        x, w = gausslegendre(T, 100)
        r = similar(x)
        dr = similar(x)
        radius_and_deriv!(r, dr, s, x)

        # Since this is actually a sphere, r should be equal everywhere, while dr should be zero.
        @test all(r .≈ 1.0)
        @test all(iszero.(dr))
    end
end
