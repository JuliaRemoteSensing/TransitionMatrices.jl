@doc raw"""
`bhmie([T=Float64,], x, m; nₘₐₓ)`

This is a simplified version of Bohren (1983)'s `bhmie` routine. Only the Mie coefficients are evaluated and returned.

Inputs:

- `T`: Type used for calculation. Real numbers will be stored as `T`, while complex numbers will be stored as `C = complex(T)`.
- `x`: Size parameter of the sphere scatterer. Defined as ``\frac{2\pi r}{\lambda}``
- `m`: Relative refractive index of the scatterer.

Keyword arguments:

- `nₘₐₓ`: Maximum order of the Mie coefficients. Default to ``\max(x + 4\sqrt{3}{x} + 2, |m|x)``.

Outputs:

- `a, b`: Mie coefficients. Both are `Vector{C}` with `nₘₐₓ` elements.

References:

- Bohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.
"""
function bhmie(T, x, m; nₘₐₓ = ceil(Int, Float64(max(x + 4 * ∛x + 2, x * abs(m)))))
    C = complex(T)
    x = T(x)
    m = C(m)
    y = m * x
    d = zeros(C, nₘₐₓ)
    for n in (nₘₐₓ - 1):-1:1
        d[n] = (n + 1) / y - (1.0 / (d[n + 1] + (n + 1) / y))
    end

    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁ = complex(ψ₁, -χ₁)
    a = zeros(C, nₘₐₓ)
    b = zeros(C, nₘₐₓ)
    for n in 1:nₘₐₓ
        ψ = (2n - 1) * ψ₁ / x - ψ₀
        χ = (2n - 1) * χ₁ / x - χ₀
        ξ = complex(ψ, -χ)
        a[n] = ((d[n] / m + n / x) * ψ - ψ₁) / ((d[n] / m + n / x) * ξ - ξ₁)
        b[n] = ((d[n] * m + n / x) * ψ - ψ₁) / ((d[n] * m + n / x) * ξ - ξ₁)
        ψ₀, ψ₁ = ψ₁, ψ
        χ₀, χ₁ = χ₁, χ
        ξ₁ = complex(ψ₁, -χ₁)
    end

    return a, b
end

bhmie(x, m; nₘₐₓ = ceil(Int, Float64(max(x + 4 * ∛x + 2, x * abs(m))))) = bhmie(Float64, x, m; nₘₐₓ = nₘₐₓ)
