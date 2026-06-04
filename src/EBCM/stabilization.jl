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
# `_F⁺_matrix` assembles the whole F⁺ matrix stably for all x (series-seeded
# recursion, §4.2); `transition_matrix(…; stable=true)` uses it to rebuild the
# spheroid `𝐔`-matrix integrands for both m=0 and m>0.

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
`(2n-1)!!`-scale cancellation. This `BigFloat` accumulation also subsumes the
additional `s≈1` cancellation that Somerville et al. 2013 handle analytically in
their Appendix B (the leading coefficients vanish at `s=1`): the base precision
has many digits of margin for it, so no special treatment is needed — the s≈1
accuracy floor of the stabilized assembly is instead the `Float64` round-off in
`_eval_F⁺`, which extra coefficient precision cannot remove.
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

The terms `c_q x^{2q+…}` first grow then decay; the sum is truncated once it has
started to decay and a term is negligible relative to the running total
(Somerville et al. 2013, §4.1). This both avoids summing past convergence and,
crucially, prevents the monomial `xp = x^{2q+…}` from overflowing to `Inf`/`NaN`
for large `x` (it would reach `x^{2·length(coeffs)}`). The series is only well
conditioned in `Float64` when the order `n` exceeds the argument `x`; in that
regime convergence is reached long before `xp` overflows.
"""
function _eval_F⁺(qmin::Integer, coeffs::AbstractVector{<:Complex}, n::Integer,
        k::Integer, x::T) where {T <: Real}
    p0 = 2qmin + k - n + 2
    x2 = x * x
    xp = x^p0
    val = zero(Complex{T})
    tol = eps(T)
    prevmag = T(Inf)
    decaying = false
    @inbounds for idx in eachindex(coeffs)
        isfinite(xp) || break
        term = Complex{T}(coeffs[idx]) * xp
        val += term
        tmag = abs(term)
        decaying |= tmag < prevmag
        if decaying && tmag ≤ tol * abs(val)
            break
        end
        prevmag = tmag
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
function _xχψ′⁺(n, k, s, x; kw...)
    ((k + 1) * F⁺(n, k - 1, s, x; kw...) - k * F⁺(n, k + 1, s, x; kw...)) / (2k + 1)
end

"""
    _xχ′ψ⁺(n, k, s, x) -> Complex

`[x·χ′_n(x)·ψ_k(sx)]⁺` (Somerville 2013, Eq. 60), from
`(2n+1)χ′_n(x) = (n+1)χ_{n-1}(x) − n χ_{n+1}(x)`.
"""
function _xχ′ψ⁺(n, k, s, x; kw...)
    ((n + 1) * F⁺(n - 1, k, s, x; kw...) - n * F⁺(n + 1, k, s, x; kw...)) / (2n + 1)
end

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

"""
    _F⁺_lastrow(s, N, xmax; prec, nterms) -> Vector

Precompute the *x-independent* power-series coefficients of the F⁺-matrix last row
`n = N+1` (the entries with `n-k ≥ 4`, `n+k` even, that seed the recursion). These
are the only expensive (`BigFloat`) part of `_F⁺_matrix`, so computing them once
and passing them to every per-quadrature-point `_F⁺_matrix` call removes the
dominant cost of the stabilized assembly. `nterms` is sized for the *largest*
argument `xmax = k·rₘₐₓ` that will be evaluated (`_eval_F⁺` truncates adaptively
for smaller `x`). Element `k+1` holds `(qmin, coeffs)`, or `nothing` where unused.
"""
function _F⁺_lastrow(s::Number, N::Integer, xmax::Real;
        prec::Integer = max(256, 12 * (N + 1) + 128),
        nterms::Integer = max(48, ceil(Int, 2 * abs(s) * xmax + 40)))
    M = N + 1
    lastrow = Vector{Union{Nothing, Tuple{Int, Vector{Complex{BigFloat}}}}}(nothing,
        M + 1)
    for k in 0:(M - 4)
        iseven(M + k) || continue
        lastrow[k + 1] = _F⁺_coeffs(M, k, s; nterms, prec)
    end
    return lastrow
end

@doc raw"""
    _F⁺_matrix(s, x, N; prec, nterms, lastrow) -> Matrix{Complex{T}}

Stably evaluate the whole matrix of `F⁺_{nk}(s,x)` for `0 ≤ n,k ≤ N+1`, returned
with a 1-based offset (`F[n+1, k+1]` holds `F⁺_{nk}`). Only `n+k`-even entries are
filled; the rest stay zero (they are only needed for `m>0`).

The bare power series is well conditioned in `Float64` *only* where the order
exceeds the argument, so it cannot be used alone for the small-`n`/large-`x`
corner. This routine combines the three regimes of Somerville, Auguié & Le Ru,
JQSRT 123 (2013) §4.2, so the result is accurate for all `x`:

  * `n ≤ k+2` — `F⁺ = F = x·χ_n(x)·ψ_k(sx)`, the direct Riccati–Bessel product
    (no negative powers ⇒ no cancellation, stable for all `x`);
  * last row `n = N+1`, deep entries `n-k ≥ 4` — the power series (Eq. 46), which
    is accurate because `N` is taken large enough that `N+1 ≫ x`;
  * `n ≥ k+4` — the stable recursion (Eq. 51, scheme (c) of their Fig. 3), filled
    along diagonals `j = n-k = 4, 6, …` from high `n` down, seeded by the last row
    and the `n = k+2` sub-diagonal:

    ```math
    F^+_{n,k} = \frac{2k+3}{s(2n+1)}\,(F^+_{n+1,k+1} + F^+_{n-1,k+1}) - F^+_{n,k+2}.
    ```

`N` must satisfy roughly `N+1 ≳ x + 15` for the last-row seed (and hence the whole
bottom-left block) to be accurate; this is comparable to the multipole order
needed for convergence anyway.
"""
function _F⁺_matrix(s::Number, x::T, N::Integer;
        prec::Integer = max(256, 12 * (N + 1) + 128),
        nterms::Integer = max(48, ceil(Int, 2 * abs(s) * x + 40)),
        lastrow = nothing) where {T <: Real}
    CT = Complex{T}
    M = N + 1
    # x-independent last-row series coefficients (precompute once across points).
    lastrow === nothing && (lastrow = _F⁺_lastrow(s, N, x; prec, nterms))
    F = zeros(CT, M + 1, M + 1)          # F[n+1, k+1] = F⁺_{nk},  n,k ∈ 0:M
    at(n, k) = @inbounds F[n + 1, k + 1]

    # Riccati–Bessel functions at this point, indexed 0:M.
    sx = s * x
    χ1, _ = ricattibessely(M, x)
    nex = estimate_ricattibesselj_extra_terms(M, abs(s) * x)
    ψ1, _ = ricattibesselj(M, nex, sx)
    χv = Vector{T}(undef, M + 1)
    ψv = Vector{CT}(undef, M + 1)
    χv[1] = -cos(x)                      # χ_0(x) = x·y_0(x) = -cos x
    ψv[1] = sin(sx)                      # ψ_0(z) = z·j_0(z) = sin z
    @inbounds for n in 1:M
        χv[n + 1] = χ1[n]
        ψv[n + 1] = ψ1[n]
    end

    # (1) direct products for n ≤ k+2 (covers all entries with no cancellation).
    @inbounds for k in 0:M, n in 0:min(k + 2, M)

        iseven(n + k) || continue
        F[n + 1, k + 1] = x * χv[n + 1] * ψv[k + 1]
    end

    # (2) last row n = M (= N+1), only the deep entries n-k ≥ 4 (the shallow part
    #     of the last row is already filled by the direct products above).
    @inbounds for k in 0:(M - 4)
        iseven(M + k) || continue
        qmin, coeffs = lastrow[k + 1]
        F[M + 1, k + 1] = _eval_F⁺(qmin, coeffs, M, k, x)
    end

    # (3) stable recursion (Eq. 51, scheme c): diagonals j = 4, 6, … filled from
    #     high n down; each entry uses the same-diagonal entry at n+1 and the j-2
    #     diagonal at (n-1,k+1) and (n,k+2), all already known.
    @inbounds for j in 4:2:M
        for n in (M - 1):-1:j
            k = n - j
            F[n + 1, k + 1] = (2k + 3) / (s * (2n + 1)) *
                              (at(n + 1, k + 1) + at(n - 1, k + 1)) - at(n, k + 2)
        end
    end

    return F
end

# ── Modified products read from a precomputed `_F⁺_matrix` ────────────────────
# The bare `_xχψ′⁺`/`_L⁷⁺`/… above call `F⁺` (the small-x-only series) directly
# and are kept for testing. Production assembly uses these matrix-backed variants
# instead, fed by `_F⁺_matrix`, which is stable for all x. `F` is the offset
# matrix returned by `_F⁺_matrix` (`F[n+1,k+1] = F⁺_{nk}`).
@inline _Fp(F, n, k) = @inbounds F[n + 1, k + 1]

"`[x·χ_n·ψ′_k]⁺` (Somerville et al., JQSRT 123 (2013), Eq. 59) from a precomputed F⁺ matrix."
_xχψ′⁺_mat(F, n, k) = ((k + 1) * _Fp(F, n, k - 1) - k * _Fp(F, n, k + 1)) / (2k + 1)

"`[x·χ′_n·ψ_k]⁺` (Somerville et al., JQSRT 123 (2013), Eq. 60) from a precomputed F⁺ matrix."
_xχ′ψ⁺_mat(F, n, k) = ((n + 1) * _Fp(F, n - 1, k) - n * _Fp(F, n + 1, k)) / (2n + 1)

"L⁷ radial factor (Somerville et al., JQSRT 123 (2013), Eq. 62) from a precomputed F⁺ matrix."
function _L⁷⁺_mat(F, n, k)
    Fmm = _Fp(F, n - 1, k - 1);
    Fpp = _Fp(F, n + 1, k + 1)
    Fmp = _Fp(F, n - 1, k + 1);
    Fpm = _Fp(F, n + 1, k - 1)
    return ((n + k + 1) * ((n + 1) * Fmm + n * Fpp) +
            (n - k) * ((n + 1) * Fmp + n * Fpm)) / ((2n + 1) * (2k + 1))
end

"L⁸ radial factor (Somerville et al., JQSRT 123 (2013), Eq. 61) from a precomputed F⁺ matrix."
function _L⁸⁺_mat(F, n, k)
    Fmm = _Fp(F, n - 1, k - 1);
    Fpp = _Fp(F, n + 1, k + 1)
    Fmp = _Fp(F, n - 1, k + 1);
    Fpm = _Fp(F, n + 1, k - 1)
    return ((n + k + 1) * ((k + 1) * Fmm + k * Fpp) +
            (k - n) * ((k + 1) * Fpm + k * Fmp)) / ((2n + 1) * (2k + 1))
end

@testitem "_F⁺_matrix is stable across x (series-seeded Eq.51 recursion)" begin
    using TransitionMatrices: _F⁺_matrix, _F⁺_coeffs, _eval_F⁺
    s = 1.5 + 0.02im
    # The bare series fails for small n / large x; the matrix builder must stay
    # accurate over the whole bottom-left (cancellation) block n ≥ k+4 for every x.
    for x in (0.5, 8.0, 16.0, 35.0)
        N = max(20, ceil(Int, x) + 18)          # N+1 ≳ x+15 for the last-row seed
        F = _F⁺_matrix(s, x, N)
        nt = max(60, ceil(Int, 2 * abs(s) * x + 50))
        maxrel = 0.0
        for n in 4:N, k in 0:(n - 4)

            iseven(n + k) || continue
            vbig = setprecision(BigFloat, 600) do
                qmin, co = _F⁺_coeffs(n, k, s; nterms = nt, prec = 600)
                _eval_F⁺(qmin, co, n, k, BigFloat(x))
            end
            ab = abs(ComplexF64(vbig))
            ab < 1e-285 && continue
            maxrel = max(maxrel,
                abs(ComplexF64(F[n + 1, k + 1]) - ComplexF64(vbig)) / ab)
        end
        @test maxrel < 1e-9
    end
end

@testitem "F⁺ matches the Bessel product where no negative powers exist (n ≤ k+2)" begin
    using TransitionMatrices: F⁺, ricattibesselj, ricattibessely,
                              estimate_ricattibesselj_extra_terms
    setprecision(BigFloat, 320) do
        s = Complex{BigFloat}(1.5, 0.02)
        for x in BigFloat.((0.3, 1.0, 2.5)),
            (n, k) in ((1, 1), (2, 2), (1, 3), (3, 5),
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
        for x in BigFloat.((0.5, 1.5)),
            (n, k) in ((1, 1), (2, 2), (3, 3), (1, 3),
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

@testitem "ebcm_matrices_m₀(; stable=true) fixes high-aspect spheroid 𝐔 vs BigFloat" begin
    using TransitionMatrices: Spheroid, ebcm_matrices_m₀, 𝐓_from_𝐏_and_𝐔

    relmax(A, B) = begin
        m = 0.0
        for I in eachindex(B)
            b = ComplexF64(B[I])
            abs(b) < 1e-12 && continue
            m = max(m, abs(ComplexF64(A[I]) - b) / abs(b))
        end
        m
    end

    λ = 2π
    a, c, N, Ng = 0.3, 5.0, 16, 300       # prolate, aspect ≈ 17: std Float64 𝐔 is hopeless
    mr = 1.5 + 0.01im
    sb = Spheroid{BigFloat, Complex{BigFloat}}(big(a), big(c), Complex{BigFloat}(mr))
    s64 = Spheroid{Float64, ComplexF64}(a, c, ComplexF64(mr))

    Pb, Ub = setprecision(BigFloat, 320) do
        ebcm_matrices_m₀(sb, big(λ), N, Ng)
    end
    Ps, Us = ebcm_matrices_m₀(s64, λ, N, Ng)
    Pst, Ust = ebcm_matrices_m₀(s64, λ, N, Ng; stable = true)

    @test relmax(Us, Ub) > 1e1            # standard products are catastrophically wrong
    @test relmax(Ust, Ub) < 1e-6          # the F⁺ matrix recovers full precision
    @test Pst == Ps                       # `stable` must leave 𝐏 untouched

    Tb = 𝐓_from_𝐏_and_𝐔(ComplexF64.(Pb), ComplexF64.(Ub))
    @test relmax(𝐓_from_𝐏_and_𝐔(Pst, Ust), Tb) < 1e-5
end

@testitem "stable=true reproduces Somerville 2013 Table 2 (prolate Qsca/Qext)" begin
    using TransitionMatrices: Spheroid, transition_matrix, scattering_cross_section,
                              extinction_cross_section

    # Somerville, Auguié & Le Ru, JQSRT 123 (2013), Table 2 (Model 1 of their
    # Ref. [35]): prolate, s = 1.55+0.01i, aspect h = 4, xmax = k·c = 10.079368.
    # Reference (their New DP / AP): Qsca = 3.21290554203156, Qext = 3.36721292620922.
    c = 10.079368
    a = c / 4
    s = Spheroid{Float64, ComplexF64}(a, c, 1.55 + 0.01im)

    # efficiencies are normalised by the equal-surface-area radius (prolate, a<c)
    e = √(1 - a^2 / c^2)
    rsurf = √(2 * a^2 * (1 + (c / (a * e)) * asin(e)) / 4)
    area = π * rsurf^2

    𝐓 = transition_matrix(s, 2π, 32, 400; stable = true)
    Qsca = scattering_cross_section(𝐓, 2π) / area
    Qext = extinction_cross_section(𝐓, 2π) / area

    @test isapprox(Qsca, 3.21290554203156; rtol = 1e-6)
    @test isapprox(Qext, 3.36721292620922; rtol = 1e-6)
end
