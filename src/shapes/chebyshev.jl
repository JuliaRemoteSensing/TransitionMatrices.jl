@doc raw"""
A Chebyshev scatterer defined by

```math
r(\theta, \phi)=r_0(1+\varepsilon T_n(\cos\theta))
```

where ``T_n(\cos\theta)=\cos n\theta``.

Attributes:

- `r₀`: Radius of the base sphere.
- `ε`: Deformation parameter, which satisfies ``-1\le\varepsilon<1``.
- `n`: Degree of the Chebyshev polynomial.
- `m`: Relative complex refractive index.
"""
struct Chebyshev{T, CT} <: AbstractAxisymmetricShape{T, CT}
    r₀::T
    ε::T
    n::Int
    m::CT
end

function volume_equivalent_radius(c::Chebyshev)
    ε = c.ε
    n = c.n
    a = 3 // 2 * ε^2 * (4n^2 - 2) / (4n^2 - 1) + 1
    if iseven(n)
        a -= 3ε * (1 + 1 // 4 * ε^2) / (n^2 - 1) + 1 // 4 * ε^3 / (9n^2 - 1)
    end

    return ∛a * c.r₀
end

volume(c::Chebyshev) = 4 // 3 * π * volume_equivalent_radius(c)^3
has_symmetric_plane(c::Chebyshev) = iseven(c.n)

@testitem "Utility functions are correct" begin
    using TransitionMatrices: Chebyshev, volume, volume_equivalent_radius,
                              has_symmetric_plane

    @testset "r₀ = $r₀, ε = $ε, n = $n" for (r₀, ε, n, V, rᵥ) in [
        (1.0, 0.5, 2, 3.4258319889145845, 0.9351741286290055),
        (3.0, -1.0, 3, 277.8963101575429, 4.048225756255873),
        (5.0, 1.0, 8, 1274.5227007709325, 6.725939997427172),
    ]
        c = Chebyshev(r₀, ε, n, 1.5 + 0.01im)
        @test volume(c) ≈ V
        @test volume_equivalent_radius(c) ≈ rᵥ
        @test has_symmetric_plane(c) == iseven(n)
    end
end

"""
```
radius_and_deriv!(r, dr, c::Chebyshev{T}, x) where {T}
```

We could have used @turbo for primitive types like `Float64`, but this is
not the bottleneck, so we only use `@inbounds` here.
"""
function radius_and_deriv!(r, dr, c::Chebyshev{T}, x) where {T}
    ngauss = length(x)

    @inbounds for i in 1:ngauss
        nϑ = acos(x[i]) * c.n
        r[i] = c.r₀ * (1 + c.ε * cos(nϑ))
        dr[i] = -c.r₀ * c.ε * c.n * sin(nϑ)
    end
end
