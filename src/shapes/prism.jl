"""
```
Prism{N, T, CT}(a, h, m)
```

A prism scatterer. The number of base edges is represented by the type parameter `N`. The prism is assumed to be aligned with the z-axis, and one of The lateral edges is assumed to have `φ=0`.


Attributes:

- `a`: the length of the prism base.
- `h`: the height of the prism.
- `m`: the relative complex refractive index.
"""
struct Prism{N, T, CT} <: AbstractNFoldShape{N, T, CT}
    a::T
    h::T
    m::CT
end

Prism(N, a, h, m) = Prism{N, typeof(a), typeof(m)}(a, h, m)

volume(p::Prism{N}) where {N} = (N * p.a^2) / (4 * tan(π / N)) * p.h
volume_equivalent_radius(p::Prism{N}) where {N} = ∛(3 * volume(p) / 4π)
has_symmetric_plane(::Prism) = true

@testitem "Utility functions are correct" begin
    using TransitionMatrices: Prism, volume, volume_equivalent_radius,
                              has_symmetric_plane

    @testset "N = $N, a = $a, h = $h" for (N, a, h, V, rᵥ) in [
        (5, 2.0, 3.0, 20.645728807067606, 1.701820965727639),
        (3, 0.5, 5.0, 0.5412658773652743, 0.5055615227215154),
        (20, 2.0, 0.1, 12.627503029350088, 1.4445845221839735),
    ]
        p = Prism(N, a, h, 1.311)
        @test volume(p) ≈ V
        @test volume_equivalent_radius(p) ≈ rᵥ
        @test has_symmetric_plane(p)
    end
end

function rmax(p::Prism{N}) where {N}
    r = p.a / (2 * sin(π / N))
    return hypot(r, p.h / 2)
end

function rmin(p::Prism{N}) where {N}
    b = p.a / (2 * tan(π / N))
    return min(b, p.h / 2)
end

function Base.:∈(x, p::Prism{4})
    return abs(x[1]) ≤ p.a / 2 && abs(x[2]) ≤ p.a / 2 && abs(x[3]) ≤ p.h / 2
end

function Base.:∈(x, p::Prism{N}) where {N}
    φ = mod2pi(atan(x[2], x[1]))
    φ′ = mod(φ, 2π / N)
    x′, y′ = Angle2d(φ′ - φ) * [x[1], x[2]]

    r′ = hypot(x′, y′)
    b = p.a / (2 * tan(π / N))
    r = b / cos(abs(φ′ - π / N))

    return r′ ≤ r && abs(x[3]) ≤ p.h / 2
end

refractive_index(p::Prism, x) = x ∈ p ? p.m : one(p.m)
