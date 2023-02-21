@doc raw"""
A cylindrical scatterer.

Attributes:

- `r`: radius of the cylinder base
- `h`: height of the cylinder
- `m`: relative complex refractive index
"""
struct Cylinder{T, CT} <: AbstractAxisymmetricShape{T, CT}
    r::T
    h::T
    m::CT
end

volume(c::Cylinder) = π * c.r^2 * c.h
volume_equivalent_radius(c::Cylinder) = ∛(3 * c.r^2 * c.h / 4)
has_symmetric_plane(::Cylinder) = true

@testitem "Utility functions are correct" begin
    using TransitionMatrices: Cylinder, volume, volume_equivalent_radius,
                              has_symmetric_plane

    @testset "r = $r, h = $h" for (r, h, V, rᵥ) in [
        (2.0, 0.5, 6.283185307179586, 1.1447142425533319),
        (1.0, 1.0, 3.141592653589793, 0.9085602964160698),
        (0.5, 2.0, 1.5707963267948966, 0.7211247851537042),
    ]
        c = Cylinder(r, h, 1.311)
        @test volume(c) ≈ V
        @test volume_equivalent_radius(c) ≈ rᵥ
        @test has_symmetric_plane(c)
    end
end

function gaussquad(c::Cylinder{T}, ngauss) where {T}
    ng = ngauss ÷ 2
    ng1 = ng ÷ 2
    ng2 = ng - ng1

    x = zeros(T, ngauss)
    w = zeros(T, ngauss)
    x1, w1 = gausslegendre(T, ng1)
    x2, w2 = gausslegendre(T, ng2)
    xx = -cos(atan(2c.r / c.h))
    @. x[1:ng1] = 0.5(xx + 1.0) * x1 + 0.5(xx - 1.0)
    @. w[1:ng1] = 0.5(xx + 1.0) * w1
    @. x[(ng1 + 1):ng] = -0.5xx * x2 + 0.5xx
    @. w[(ng1 + 1):ng] = -0.5xx * w2
    @. x[(ng + 1):ngauss] = (-1.0) * x[ng:-1:1]
    @. w[(ng + 1):ngauss] = w[ng:-1:1]

    h = c.h / 2
    d = c.r
    r = similar(x)
    r′ = similar(x)

    # We could have used @turbo for primitive types like 
    # `Float64`, but this is not the bottleneck, so we only use 
    # `@simd` here.
    @simd for i in 1:(ngauss ÷ 2)
        cosϑ = abs(x[i])
        sinϑ = √(1 - cosϑ^2)
        if h / cosϑ < d / sinϑ
            r[i] = h / cosϑ
            r′[i] = h * sinϑ / cosϑ^2
        else
            r[i] = d / sinϑ
            r′[i] = -d * cosϑ / sinϑ^2
        end
        r[ngauss + 1 - i] = r[i]
        r′[ngauss + 1 - i] = r′[i]
        r′[i] = -r′[i]
    end

    return x, w, r, r′
end
