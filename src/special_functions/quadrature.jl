@doc raw"""
`gausslegendre([T=Float64,], n)`

Calculate the `n`-point Gauss-Legendre quadrature nodes and weights.

- For `Float64`, we use the ``\mathcal{O}(n)`` implementation from `FastGaussQuadrature.jl`.
- For other types, we use the arbitrary precision implementation from `Arblib.jl`, then convert the results to the desired type.
"""
gausslegendre(::Type{Float64}, n::Integer) = FastGaussQuadrature.gausslegendre(n)
gausslegendre(n::Integer) = gausslegendre(Float64, n)

function gausslegendre(T, n::Integer)
    x, w = gausslegendre(Arb, n)
    x = convert.(T, BigFloat.(x))
    w = convert.(T, BigFloat.(w))
    return x, w
end

function gausslegendre(::Type{<:ArbLike}, n::Integer)
    x = ArbRefVector(n)
    w = ArbRefVector(n)
    gausslegendre!(x, w, n)
    return x, w
end

function gausslegendre!(x, w, n::Integer)
    for i in 1:(n ÷ 2)
        Arblib.hypgeom_legendre_p_ui_root!(x[n + 1 - i], w[n + 1 - i], UInt64(n),
                                           UInt64(i - 1))
        x[i] = -x[n + 1 - i]
        w[i] = w[n + 1 - i]
    end

    if n % 2 == 1
        Arblib.hypgeom_legendre_p_ui_root!(x[n ÷ 2 + 1], w[n ÷ 2 + 1], UInt64(n),
                                           UInt64(n ÷ 2))
    end
end

@testitem "hypgeom_legendre_p_ui_root! is used correctly" begin
    using Arblib: Arb
    using TransitionMatrices: gausslegendre

    @testset "n = $n" for n in (2, 5, 10, 20, 50, 100, 101)
        x, w = gausslegendre(n)
        x₁, w₁ = gausslegendre(Arb, n)

        @test all(x .≈ x₁)
        @test all(w .≈ w₁)
    end
end
