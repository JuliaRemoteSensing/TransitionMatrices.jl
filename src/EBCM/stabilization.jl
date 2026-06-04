# ── Numerically stable EBCM integrands for spheroids (Somerville/Le Ru) ──────
#
# The EBCM `U`-matrix integrands contain products χ_n(x)·ψ_k(sx) of an
# irregular Riccati–Bessel function χ_n (which behaves like x^{-n} for small x)
# and a regular one ψ_k. For a spheroid, the negative-power part of their
# Laurent expansion integrates *exactly* to zero, but evaluated directly it
# dominates the integrand by many orders of magnitude and destroys precision in
# `Float64` (Somerville, Auguié & Le Ru, JQSRT 113:524 (2012)).
#
# The fix (Somerville et al., JQSRT 123:153 (2013)) replaces χ_n(x)·ψ_k(sx) by
# F⁺_{nk}(s,x)/x, where (their Eq. 45)
#
#     F_{nk}(s,x)  = x · χ_n(x) · ψ_k(sx),
#     F⁺_{nk}(s,x) = P⁺[F_{nk}]   (the part with non-negative powers of x),
#
# computed directly from the power series so the cancelling negative powers are
# never formed (their Eq. 46 and §4.1).
#
# Power series (DLMF 10.53.1/10.53.2). This package uses ψ_k = z·j_k and the
# convention χ_n = z·y_n (NO minus sign; see `ricattibessely`):
#
#     ψ_k(z) = Σ_{b≥0} β_{k,b} z^{k+1+2b},  β_{k,b} = (-1)^b / (2^b b! (2k+2b+1)!!)
#     χ_n(z) = Σ_{a≥0} γ_{n,a} z^{2a-n},
#         γ_{n,a} = -(2n-2a-1)!! / (2^a a!)                      (a ≤ n)
#         γ_{n,a} = (-1)^{n+a+1} / (2^a a! (2a-2n-1)!!)          (a ≥ n+1)
#
# Hence F_{nk}(s,x) = Σ_{q≥0} c_{nk,q}(s) x^{2q+k-n+2}, with
#
#     c_{nk,q}(s) = Σ_{a=0}^{q} γ_{n,a} β_{k,q-a} s^{k+1+2(q-a)},
#
# and F⁺ keeps q ≥ qmin = max(0, (n-k-2)÷2) (the blocks needed have n+k even, so
# n-k-2 is even). The leading coefficients c_{nk,q} themselves suffer
# catastrophic cancellation (γ_{n,0} = (2n-1)!! is enormous and the β alternate
# in sign), so they are accumulated in high-precision `BigFloat` once; the
# resulting F⁺ is a well-behaved series that is cheap and safe to evaluate in
# the requested floating-point type.
#
# This module is self-contained and not yet wired into the EBCM solver; it is
# the foundation (stage S1) for a stabilised spheroid integrand.

"""
    _odd_dfact(m) -> BigInt

Odd double factorial `m!!` for odd `m ≥ -1`, with the convention `(-1)!! = 1`.
"""
function _odd_dfact(m::Integer)
    isodd(m) || throw(ArgumentError("_odd_dfact expects an odd integer, got $m"))
    m ≥ -1 || throw(ArgumentError("_odd_dfact expects m ≥ -1, got $m"))
    r = BigInt(1)
    j = BigInt(m)
    while j > 1
        r *= j
        j -= 2
    end
    return r
end

"""
    _γ(n, a) -> Rational{BigInt}

`γ_{n,a}`: coefficient of `z^{2a-n}` in the series of `χ_n(z) = -z y_n(z)`
(DLMF 10.53.2).
"""
function _γ(n::Integer, a::Integer)
    # This package uses χ_n(z) = z·y_n(z) (NO minus sign; see `ricattibessely`),
    # which is the negative of the DLMF/Mishchenko χ = -z·y_n. Hence
    #   χ_n(z) = (-1)^{n+1} Σ_a (-1)^a z^{2a-n} / (2^a a! (2a-2n-1)!!),
    # with the negative odd double-factorial convention folded into the two
    # branches below (a ≤ n: singular/principal part, a ≥ n+1: regular tail).
    den = (BigInt(2)^a) * Base.factorial(big(a))   # Base: the package overrides `factorial`
    if a ≤ n
        return -_odd_dfact(2n - 2a - 1) // den
    else
        sgn = iseven(n + a) ? big(-1) : big(1)
        return sgn // (den * _odd_dfact(2a - 2n - 1))
    end
end

"""
    _β(k, b) -> Rational{BigInt}

`β_{k,b}`: coefficient of `z^{k+1+2b}` in the series of `ψ_k(z) = z j_k(z)`
(DLMF 10.53.1).
"""
function _β(k::Integer, b::Integer)
    sgn = iseven(b) ? big(1) : big(-1)
    den = (BigInt(2)^b) * Base.factorial(big(b)) * _odd_dfact(2k + 2b + 1)
    return sgn // den
end

"""
    _F⁺_coeffs(n, k, s; nterms, prec) -> (qmin, coeffs)

High-precision coefficients of `F⁺_{nk}(s,x) = Σ c_q x^{2q+k-n+2}` for
`q = qmin, …, qmin+nterms-1`, returned as a `Vector{Complex{BigFloat}}`.
`s` is the (generally complex) relative refractive index. The sum defining each
`c_q` is accumulated in `BigFloat` at precision `prec` bits to defeat the
`(2n-1)!!`-scale cancellation.
"""
function _F⁺_coeffs(n::Integer, k::Integer, s::Number; nterms::Integer = 48,
                       prec::Integer = max(256, 12n + 128))
    qmin = max(0, (n - k - 2) ÷ 2)
    qmax = qmin + nterms - 1
    coeffs = Vector{Complex{BigFloat}}(undef, nterms)
    setprecision(BigFloat, prec) do
        sb = Complex{BigFloat}(s)
        γ = [BigFloat(_γ(n, a)) for a in 0:qmax]
        β = [BigFloat(_β(k, b)) for b in 0:qmax]
        spow = [sb^(k + 1 + 2j) for j in 0:qmax]   # s^{k+1+2(q-a)}, indexed by (q-a)
        for (idx, q) in enumerate(qmin:qmax)
            acc = zero(Complex{BigFloat})
            for a in 0:q
                acc += γ[a + 1] * β[q - a + 1] * spow[q - a + 1]
            end
            coeffs[idx] = acc
        end
    end
    return qmin, coeffs
end

"""
    _eval_F⁺(qmin, coeffs, n, k, x::T) -> Complex{T}

Evaluate `F⁺_{nk}(s,x) = Σ_q c_q x^{2q+k-n+2}` at `x` in floating type `T`,
converting the high-precision `coeffs` to `Complex{T}` first.
"""
function _eval_F⁺(qmin::Integer, coeffs::AbstractVector{<:Complex}, n::Integer,
                     k::Integer, x::T) where {T <: Real}
    p0 = 2qmin + k - n + 2
    x2 = x * x
    xp = x^p0
    val = zero(Complex{T})
    @inbounds for idx in eachindex(coeffs)
        val += Complex{T}(coeffs[idx]) * xp
        xp *= x2
    end
    return val
end

"""
    F⁺(n, k, s, x; nterms, prec) -> Complex

Numerically stable `F⁺_{nk}(s,x)` (Somerville et al. 2013, Eq. 45–46),
evaluated in the precision of `x`. `F⁺/x` is the cancellation-free replacement
for `χ_n(x)·ψ_k(sx)` in the spheroid EBCM `U`-matrix integrand.
"""
function F⁺(n::Integer, k::Integer, s::Number, x::T;
               nterms::Integer = max(24, ceil(Int, 2 * abs(s) * x + 16)),
               prec::Integer = max(256, 12n + 128)) where {T <: Real}
    qmin, coeffs = _F⁺_coeffs(n, k, s; nterms, prec)
    return _eval_F⁺(qmin, coeffs, n, k, x)
end

# ── Cancellation-free "modified" Bessel products via F⁺ (Somerville 2013) ─────
# Each returns the P⁺ (non-negative-power) part of a product appearing in the
# spheroid EBCM U-matrix integrand, assembled from F⁺ at shifted indices through
# the Riccati–Bessel recurrences (their Eqs. 57–62). All carry the leading `x`
# factor, like `F⁺` itself. For small enough n−k these equal the full products
# (no negative powers), which the tests exploit for validation.

"""
    _xχψ′⁺(n, k, s, x) -> Complex

`[x·χ_n(x)·ψ′_k(sx)]⁺` (Somerville 2013, Eq. 59), from
`(2k+1)ψ′_k(sx) = (k+1)ψ_{k-1}(sx) − k ψ_{k+1}(sx)`.
"""
_xχψ′⁺(n, k, s, x; kw...) =
    ((k + 1) * F⁺(n, k - 1, s, x; kw...) - k * F⁺(n, k + 1, s, x; kw...)) / (2k + 1)

"""
    _xχ′ψ⁺(n, k, s, x) -> Complex

`[x·χ′_n(x)·ψ_k(sx)]⁺` (Somerville 2013, Eq. 60), from
`(2n+1)χ′_n(x) = (n+1)χ_{n-1}(x) − n χ_{n+1}(x)`.
"""
_xχ′ψ⁺(n, k, s, x; kw...) =
    ((n + 1) * F⁺(n - 1, k, s, x; kw...) - n * F⁺(n + 1, k, s, x; kw...)) / (2n + 1)

"""
    _L⁷⁺(n, k, s, x) -> Complex

`[x·(χ′_n ψ′_k + n(n+1) χ_n ψ_k /(s x²))]⁺` (Somerville 2013, Eq. 62) — the
cancellation-free radial factor of the `L⁷` integrand.
"""
function _L⁷⁺(n, k, s, x; kw...)
    Fmm = F⁺(n - 1, k - 1, s, x; kw...)
    Fpp = F⁺(n + 1, k + 1, s, x; kw...)
    Fmp = F⁺(n - 1, k + 1, s, x; kw...)
    Fpm = F⁺(n + 1, k - 1, s, x; kw...)
    return ((n + k + 1) * ((n + 1) * Fmm + n * Fpp) +
            (n - k) * ((n + 1) * Fmp + n * Fpm)) / ((2n + 1) * (2k + 1))
end

"""
    _L⁸⁺(n, k, s, x) -> Complex

`[x·(χ′_n ψ′_k + k(k+1) χ_n ψ_k /(s x²))]⁺` (Somerville 2013, Eq. 61) — the
cancellation-free radial factor of the `L⁸` integrand.
"""
function _L⁸⁺(n, k, s, x; kw...)
    Fmm = F⁺(n - 1, k - 1, s, x; kw...)
    Fpp = F⁺(n + 1, k + 1, s, x; kw...)
    Fmp = F⁺(n - 1, k + 1, s, x; kw...)
    Fpm = F⁺(n + 1, k - 1, s, x; kw...)
    return ((n + k + 1) * ((k + 1) * Fmm + k * Fpp) +
            (k - n) * ((k + 1) * Fpm + k * Fmp)) / ((2n + 1) * (2k + 1))
end

@testitem "F⁺ matches the Bessel product where no negative powers exist (n ≤ k+2)" begin
    using TransitionMatrices: F⁺, ricattibesselj, ricattibessely,
                              estimate_ricattibesselj_extra_terms
    setprecision(BigFloat, 320) do
        s = Complex{BigFloat}(1.5, 0.02)
        for x in BigFloat.((0.3, 1.0, 2.5)), (n, k) in ((1, 1), (2, 2), (1, 3), (3, 5),
                                                        (4, 4), (2, 4))
            # n ≤ k+2  ⇒  F⁺ = F = x·χ_n(x)·ψ_k(s·x) (the full product, no truncation)
            @assert n ≤ k + 2
            # request ≥ 2 orders: ricattibessely! unconditionally writes χ[2].
            χ, _ = ricattibessely(max(n, 2), x)
            nextra = estimate_ricattibesselj_extra_terms(max(k, 2), abs(s) * x)
            ψ, _ = ricattibesselj(max(k, 2), nextra, s * x)
            ref = x * χ[n] * ψ[k]
            got = F⁺(n, k, s, x; nterms = 80)
            @test isapprox(got, ref; rtol = 1e-25, atol = 1e-25)
        end
    end
end

@testitem "F⁺ satisfies the Somerville recurrence Eq. 51" begin
    using TransitionMatrices: F⁺
    # F⁺_{n+1,k} + F⁺_{n-1,k} = s (2n+1)/(2k+1) (F⁺_{n,k+1} + F⁺_{n,k-1})
    setprecision(BigFloat, 512) do
        s = Complex{BigFloat}(1.5, 0.02)
        for x in BigFloat.((0.1, 0.5, 2.0)), (n, k) in ((8, 3), (10, 2), (12, 6), (15, 4))
            nt = 80
            lhs = F⁺(n + 1, k, s, x; nterms = nt) + F⁺(n - 1, k, s, x; nterms = nt)
            rhs = s * (2n + 1) / (2k + 1) *
                  (F⁺(n, k + 1, s, x; nterms = nt) + F⁺(n, k - 1, s, x; nterms = nt))
            @test isapprox(lhs, rhs; rtol = 1e-30, atol = 1e-30)
        end
    end
end

@testitem "F⁺ Float64 evaluation matches BigFloat (no precision loss in the kept series)" begin
    using TransitionMatrices: _F⁺_coeffs, _eval_F⁺
    s = ComplexF64(1.5, 0.02)
    for (n, k) in ((10, 2), (16, 4), (20, 6)), x64 in (0.1, 0.5, 1.0, 2.0)
        qmin, coeffs = _F⁺_coeffs(n, k, s; nterms = 60, prec = 600)
        v64 = _eval_F⁺(qmin, coeffs, n, k, x64)
        vbig = _eval_F⁺(qmin, coeffs, n, k, BigFloat(x64))
        @test isapprox(ComplexF64(vbig), v64; rtol = 1e-12, atol = 1e-300)
    end
end

@testitem "Modified Bessel products (Eq. 59-62) match direct products for n ≤ k" begin
    using TransitionMatrices: _xχψ′⁺, _xχ′ψ⁺, _L⁷⁺, _L⁸⁺, ricattibesselj,
                              ricattibessely, estimate_ricattibesselj_extra_terms
    setprecision(BigFloat, 320) do
        s = Complex{BigFloat}(1.5, 0.02)
        for x in BigFloat.((0.5, 1.5)), (n, k) in ((1, 1), (2, 2), (3, 3), (1, 3),
                                                   (2, 4), (3, 5), (2, 5), (4, 6))
            # n ≤ k ⇒ all four products have only non-negative powers ⇒ ⁺ = full.
            @assert n ≤ k
            χ, χ′ = ricattibessely(max(n, 2), x)
            nextra = estimate_ricattibesselj_extra_terms(max(k, 2), abs(s) * x)
            ψ, ψ′ = ricattibesselj(max(k, 2), nextra, s * x)
            nt = 80
            @test isapprox(_xχψ′⁺(n, k, s, x; nterms = nt), x * χ[n] * ψ′[k];
                           rtol = 1e-22, atol = 1e-22)
            @test isapprox(_xχ′ψ⁺(n, k, s, x; nterms = nt), x * χ′[n] * ψ[k];
                           rtol = 1e-22, atol = 1e-22)
            l7 = x * (χ′[n] * ψ′[k] + n * (n + 1) * χ[n] * ψ[k] / (s * x^2))
            l8 = x * (χ′[n] * ψ′[k] + k * (k + 1) * χ[n] * ψ[k] / (s * x^2))
            @test isapprox(_L⁷⁺(n, k, s, x; nterms = nt), l7; rtol = 1e-22, atol = 1e-22)
            @test isapprox(_L⁸⁺(n, k, s, x; nterms = nt), l8; rtol = 1e-22, atol = 1e-22)
        end
    end
end
