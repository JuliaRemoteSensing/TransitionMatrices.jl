@doc raw"""
`bhcoat([T=Float64,], x_core, x_mantle, m_core, m_mantle; nₘₐₓ, tolerance = 1e-8)`

Inputs:

- `T`: Type used for calculation. All real numbers will be stored as `T`, while complex numbers will be stored as `C = complex(T)`.
- `x_core`: Size parameter of the core. Defined as ``\frac{2\pi r}{\lambda}``
- `x_mantle`: Size parameter of the coated sphere. `x_mantle >= x_core` should hold.
- `m_core`: Refractive index of the core, relative to the host medium.
- `m_mantle`: Refractive index of the mantle, relative to the host medium.

Keyword arguments:

- `nₘₐₓ`: Maximum order of the Mie coefficients. Default to ``\max(x_m + 4\sqrt{3}{x_m} + 2, \max(|m_c|, |m_m|)x_m)``.
- `tolerance`: Error tolerance. Default is `1e-8`.

Outputs:

- `a, b`: Mie coefficients. Both are of type `Vector{C}`.

References:

- Bohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.
"""
function bhcoat(T, x_core, x_mantle, m_core, m_mantle;
                nₘₐₓ = ceil(Int,
                            Float64(max(x_mantle + 4 * ∛x_mantle + 2,
                                        x_mantle * max(abs(m_core), abs(m_mantle))))),
                tolerance = 1e-8)
    @assert x_mantle>=x_core "x_mantle must be greater than or equal to x_core"

    C = complex(T)
    x = T(x_core)
    y = T(x_mantle)
    m₁ = C(m_core)
    m₂ = C(m_mantle)
    x₁ = m₁ * x
    x₂ = m₂ * x
    y₂ = m₂ * y
    m = m₂ / m₁
    a = zeros(C, nₘₐₓ)
    b = zeros(C, nₘₐₓ)
    d0x1 = cot(x₁)
    d0x2 = cot(x₂)
    d0y2 = cot(y₂)
    d1x1 = zero(T)
    d1x2 = zero(T)
    brack = zero(T)
    crack = zero(T)
    psi0y = cos(y)
    psi1y = sin(y)
    chi0y = -sin(y)
    chi1y = cos(y)
    xi0y = psi0y - 1.0im * chi0y
    xi1y = psi1y - 1.0im * chi1y
    chi0y2 = -sin(y₂)
    chi1y2 = cos(y₂)
    chipy2 = zero(T)
    chiy2 = zero(T)
    chi0x2 = -sin(x₂)
    chi1x2 = cos(x₂)
    chipx2 = zero(T)
    chix2 = zero(T)
    amess1 = zero(T)
    amess2 = zero(T)
    amess3 = zero(T)
    amess4 = zero(T)
    iflag = false
    for n in 1:nₘₐₓ
        rn = T(n)
        psiy = (2rn - 1) * psi1y / y - psi0y
        chiy = (2rn - 1) * chi1y / y - chi0y
        xiy = psiy - 1.0im * chiy
        d1y2 = 1 / (rn / y₂ - d0y2) - rn / y₂
        if !iflag
            d1x1 = 1 / (rn / x₁ - d0x1) - rn / x₁
            d1x2 = 1 / (rn / x₂ - d0x2) - rn / x₂
            chix2 = (2rn - 1) * chi1x2 / x₂ - chi0x2
            chiy2 = (2rn - 1) * chi1y2 / y₂ - chi0y2
            chipx2 = chi1x2 - rn * chix2 / x₂
            chipy2 = chi1y2 - rn * chiy2 / y₂
            ancap = m * d1x1 - d1x2
            ancap = ancap / (m * d1x1 * chix2 - chipx2)
            ancap = ancap / (chix2 * d1x2 - chipx2)
            brack = ancap * (chiy2 * d1y2 - chipy2)
            bncap = m * d1x2 - d1x1
            bncap = bncap / (m * chipx2 - d1x1 * chix2)
            bncap = bncap / (chix2 * d1x2 - chipx2)
            crack = bncap * (chiy2 * d1y2 - chipy2)
            amess1 = brack * chipy2
            amess2 = brack * chiy2
            amess3 = crack * chipy2
            amess4 = crack * chiy2
        end
        if abs(amess1) < tolerance * abs(d1y2) && abs(amess2) < tolerance &&
           abs(amess3) < tolerance * abs(d1y2) && abs(amess4) < tolerance
            brack = zero(T)
            crack = zero(T)
            iflag = true
        else
            iflag = false
        end
        dnbar = d1y2 - brack * chipy2
        dnbar = dnbar / (1.0 - brack * chiy2)
        gnbar = d1y2 - crack * chipy2
        gnbar = gnbar / (1.0 - crack * chiy2)
        a[n] = ((dnbar / m₂ + rn / y) * psiy - psi1y) / ((dnbar / m₂ + rn / y) * xiy - xi1y)
        b[n] = ((m₂ * gnbar + rn / y) * psiy - psi1y) / ((m₂ * gnbar + rn / y) * xiy - xi1y)
        psi0y = psi1y
        psi1y = psiy
        chi0y = chi1y
        chi1y = chiy
        xi1y = psi1y - 1.0im * chi1y
        chi0x2 = chi1x2
        chi1x2 = chix2
        chi0y2 = chi1y2
        chi1y2 = chiy2
        d0x1 = d1x1
        d0x2 = d1x2
        d0y2 = d1y2
    end

    return a, b
end

function bhcoat(x_core, x_mantle, m_core, m_mantle;
                nₘₐₓ = ceil(Int,
                            Float64(max(x_mantle + 4 * ∛x_mantle + 2,
                                        x_mantle * max(abs(m_core), abs(m_mantle))))),
                tolerance = 1e-8)
    bhcoat(Float64, x_core, x_mantle, m_core, m_mantle; nₘₐₓ = nₘₐₓ, tolerance = tolerance)
end

@testitem "bhcoat" begin
    using TransitionMatrices: bhmie, bhcoat

    @testset "converges to bhmie when x = $x, m = $m" for (x, m) in [
        (1.0, 1.311),
        (2.0, 1.5 + 0.01im),
    ]
        am, bm = bhmie(x, m)
        ac, bc = bhcoat(x, x, m, m)
        @test all(isapprox.(am, ac; atol = 1e-12))
        @test all(isapprox.(bm, bc; atol = 1e-12))
    end
end
