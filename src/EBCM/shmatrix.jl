# ── Sh-matrix moment-separation feature ──────────────────────────────────────
#
# The EBCM surface integrals that build `𝐏` and `𝐔` all have the form
#
#     ∫ (angular×geometry) · f_n(kr) · g_{n′}(s·kr) dϑ ,
#
# where `f,g` are Riccati–Bessel functions `ψ,ψ′,χ,χ′`, `k = 2π/λ` is the size
# parameter and `s` the relative refractive index. Expanding the radial
# functions in their power series (DLMF 10.53, coefficients `_β`,`_γ` — the same
# tables used by the `F⁺` stabilization) factors every integral as
#
#     ∫ = Σ_{terms} [size/material coefficient in (k,s)] × [shape-only moment M(q)] ,
#
# where the *shape moment* `M(q) = ∫ (angular×geometry) · r^q dϑ` depends only on
# the particle geometry and the azimuthal index `m`, **not** on `k` or `s`. The
# moments are therefore computed once (`prepare_sh`); each `(λ, s)` evaluation is
# a cheap coefficient×moment sum, which is what makes parameter sweeps fast.
#
# This is the same idea as the `F⁺` construction (expand and integrate term by
# term), generalized to every EBCM integral. It also recovers the spheroid `𝐔`
# stabilization for free: for a spheroid the negative-`r`-power moments vanish
# analytically, so — computed at high precision — they contribute ≈0 and the
# catastrophic cancellation in the irregular `χ_n·ψ_{n′}` products never forms.
#
# References:
# - V. G. Farafonov, Sh-matrix family (parameter separation).
# - Somerville, Auguié & Le Ru, JQSRT 123 (2013) 153 (the `F⁺` projection).

# ── Radial power-series term tables (reuse `_β`, `_γ`) ────────────────────────
#
# For order `n`, each Riccati–Bessel function is `Σ_{b≥0} c_b · z^{base+2b}`:
#   ψ_n(z)  : base = n+1,   c_b = β_{n,b}
#   ψ′_n(z) : base = n,     c_b = β_{n,b}·(n+1+2b)
#   χ_n(z)  : base = -n,    c_b = γ_{n,b}
#   χ′_n(z) : base = -n-1,  c_b = γ_{n,b}·(2b-n)
"""
    _sh_series(kind, n, B, R) -> (base, coeffs::Vector{R})

Coefficients `c_b` (`b = 0:B-1`) and the lowest exponent `base` of the power
series of a Riccati–Bessel function `kind ∈ (:ψ, :ψ′, :χ, :χ′)` of order `n`,
i.e. `f_n(z) = Σ_b coeffs[b+1] · z^{base + 2b}`.
"""
function _sh_series(kind::Symbol, n::Integer, B::Integer, ::Type{R}) where {R}
    coeffs = Vector{R}(undef, B)
    if kind === :ψ
        base = n + 1
        for b in 0:(B - 1)
            coeffs[b + 1] = R(_β(n, b))
        end
    elseif kind === :ψ′
        base = n
        for b in 0:(B - 1)
            coeffs[b + 1] = R(_β(n, b)) * (n + 1 + 2b)
        end
    elseif kind === :χ
        base = -n
        for b in 0:(B - 1)
            coeffs[b + 1] = R(_γ(n, b))
        end
    elseif kind === :χ′
        base = -n - 1
        for b in 0:(B - 1)
            coeffs[b + 1] = R(_γ(n, b)) * (2b - n)
        end
    else
        throw(ArgumentError("unknown Riccati–Bessel kind $kind"))
    end
    return base, coeffs
end

# ── Coefficient×moment reconstruction ─────────────────────────────────────────
"""
    _sh_recon(fterm, gterm, k, s, M, shift) -> Complex

Reconstruct `∫ (geometry) · f_n(kr) · g_{n′}(s·kr) · r^{shift} dϑ` from the
radial series terms `fterm = (basef, cf)`, `gterm = (baseg, cg)` and the
shape-moment lookup `M(q)` (returns the moment of `r^q`). The size/material
dependence is entirely in `k` (real) and `s` (complex); `M` carries the
geometry. Powers are advanced incrementally to avoid `^` in the inner loop.
"""
function _sh_recon(fterm, gterm, k::Real, s::Number, M, shift::Integer)
    basef, cf = fterm
    baseg, cg = gterm
    Bf = length(cf)
    Bg = length(cg)
    k2 = k * k
    sk = s * k
    s2k2 = sk * sk
    pref = k^basef * sk^baseg
    acc = zero(typeof(pref))
    k2b = one(k)
    @inbounds for b in 0:(Bf - 1)
        s2c = one(s2k2)
        cfb = cf[b + 1]
        for c in 0:(Bg - 1)
            q = basef + baseg + shift + 2 * (b + c)
            acc += (cfb * cg[c + 1]) * (k2b * s2c) * M(q)
            s2c *= s2k2
        end
        k2b *= k2
    end
    return pref * acc
end

# ── High-precision shape widening (for the vanishing negative-power moments) ──
function _sh_widen(s::Spheroid, ::Type{R}) where {R}
    Spheroid{R, Complex{R}}(R(s.a), R(s.c), Complex{R}(s.m))
end
# Fallback: a general shape is used at its native precision. The negative-power
# moments then carry only native round-off, so the analytic `𝐔` stabilization is
# *not* available for non-spheroid shapes (gated in `prepare_sh`); `𝐏` and the
# parameter-sweep reuse work for any axisymmetric shape.
_sh_widen(s, ::Type{R}) where {R} = s

# ── Shape moments for one azimuthal index `m` ────────────────────────────────
"""
    _sh_moments_m(shape, m, nmax, Ng, qlo, qhi; momtype, store) -> NamedTuple

Compute the shape-only moment tables for azimuthal index `m`, accumulated at
precision `momtype` (default `BigFloat`, so the spheroid's vanishing
negative-power moments are accurate) and stored as `store` (default `Float64`).
Families (sum over `i in 1:ng`, mirroring `ebcm_matrices_m₀`):

  Mτd[n,n′,q]  = Σ wᵢ r′ᵢ τ[i,n] d[i,n′] rᵢ^q
  Mdτ[n,n′,q]  = Σ wᵢ r′ᵢ d[i,n] τ[i,n′] rᵢ^q
  Mπd[n,n′,q]  = Σ wᵢ r′ᵢ π[i,n] d[i,n′] rᵢ^q   (only `m>0`)
  Mππττ[n,q]   = Σ wᵢ (π[i,n]²+τ[i,n]²) rᵢ^q     (no r′)
"""
function _sh_moments_m(shape, m::Integer, nmax::Integer, Ng::Integer, qlo::Integer,
        qhi::Integer; momtype::Type{R} = BigFloat,
        store::Type{Ts} = Float64) where {R, Ts}
    sym = has_symmetric_plane(shape)
    ng = sym ? Ng ÷ 2 : Ng
    nmin = max(1, Int(m))

    sh = _sh_widen(shape, R)
    x, w, r, r′ = gaussquad(sh, Ng)
    ϑ = acos.(x)

    d = OffsetArray(zeros(R, Ng, nmax + 1), 1:Ng, 0:nmax)
    τ = similar(d)
    𝜋 = similar(d)
    for i in 1:Ng
        wigner_d_recursion!(view(d, i, :), 0, Int(m), nmax, ϑ[i]; deriv = view(τ, i, :))
        for n in 0:nmax
            𝜋[i, n] = pi_func(R, Int(m), n, ϑ[i]; d = d[i, n])
        end
    end

    nq = qhi - qlo + 1
    rq = OffsetArray(zeros(R, ng, nq), 1:ng, qlo:qhi)
    for i in 1:ng
        ri = r[i]
        rq[i, 0] = one(R)
        for q in 1:qhi
            rq[i, q] = rq[i, q - 1] * ri
        end
        for q in -1:-1:qlo
            rq[i, q] = rq[i, q + 1] / ri
        end
    end

    Mτd = OffsetArray(zeros(Ts, nmax, nmax, nq), 1:nmax, 1:nmax, qlo:qhi)
    Mdτ = OffsetArray(zeros(Ts, nmax, nmax, nq), 1:nmax, 1:nmax, qlo:qhi)
    Mπd = OffsetArray(zeros(Ts, nmax, nmax, nq), 1:nmax, 1:nmax, qlo:qhi)
    Mππττ = OffsetArray(zeros(Ts, nmax, nq), 1:nmax, qlo:qhi)

    for n in nmin:nmax, n′ in nmin:nmax

        for q in qlo:qhi
            sτd = zero(R)
            sdτ = zero(R)
            sπd = zero(R)
            for i in 1:ng
                wr = w[i] * r′[i] * rq[i, q]
                sτd += wr * τ[i, n] * d[i, n′]
                sdτ += wr * d[i, n] * τ[i, n′]
                sπd += wr * 𝜋[i, n] * d[i, n′]
            end
            Mτd[n, n′, q] = Ts(sτd)
            Mdτ[n, n′, q] = Ts(sdτ)
            Mπd[n, n′, q] = Ts(sπd)
        end
    end
    for n in nmin:nmax, q in qlo:qhi

        sm = zero(R)
        for i in 1:ng
            sm += w[i] * (𝜋[i, n]^2 + τ[i, n]^2) * rq[i, q]
        end
        Mππττ[n, q] = Ts(sm)
    end

    return (; Mτd, Mdτ, Mπd, Mππττ, nmin, qlo, qhi, sym, ng)
end

# Default `q` band wide enough for every `m=0..nmax` integrand:
#   lowest  : χ′_n·ψ_{n′}/(kr)²  → (-n-1)+(n′)-2 ≥ -(nmax+2)
#   highest : ψ_n·ψ_{n′}         → (n+1)+(n′+1)+2(B-1) ≤ 2nmax+2B
_sh_qband(nmax::Integer, B::Integer) = (-(nmax + 2), 2nmax + 2B)

# Bounds-safe moment lookup (terms outside the stored band are negligibly small
# high-order tails / analytically absent low orders).
@inline function _sh_lookup(M3, n::Integer, n′::Integer, q::Integer, qlo::Integer,
        qhi::Integer)
    (qlo ≤ q ≤ qhi) ? (@inbounds M3[n, n′, q]) : zero(eltype(M3))
end
@inline function _sh_lookup(M2, n::Integer, q::Integer, qlo::Integer, qhi::Integer)
    (qlo ≤ q ≤ qhi) ? (@inbounds M2[n, q]) : zero(eltype(M2))
end

# ── m=0 assembly from moments (mirrors `ebcm_matrices_m₀`) ────────────────────
"""
    _sh_matrices_m₀(mom, k, s, nmax, B, CT) -> (𝐏, 𝐔)

Reconstruct the `m=0` `𝐏` and `𝐔` blocks from the precomputed `m=0` moment
tables `mom` at size parameter `k` and refractive index `s`, using `B` radial
series terms. The algebra reproduces `ebcm_matrices_m₀` term by term.
"""
function _sh_matrices_m₀(mom, k::Real, s::Number, nmax::Integer, B::Integer,
        ::Type{CT}) where {CT}
    R = real(CT)
    sym = mom.sym
    qlo, qhi = mom.qlo, mom.qhi
    Mτd, Mdτ, Mππττ = mom.Mτd, mom.Mdτ, mom.Mππττ

    a = [n * (n + 1) for n in 1:nmax]
    A = [√(R(2n + 1) / (2n * (n + 1))) for n in 1:nmax]

    ψ = [_sh_series(:ψ, n, B, R) for n in 1:nmax]
    ψ′ = [_sh_series(:ψ′, n, B, R) for n in 1:nmax]
    χ = [_sh_series(:χ, n, B, R) for n in 1:nmax]
    χ′ = [_sh_series(:χ′, n, B, R) for n in 1:nmax]

    𝐏 = zeros(CT, 2nmax, 2nmax)
    𝐏₁₁ = view(𝐏, 1:nmax, 1:nmax)
    𝐏₂₂ = view(𝐏, (nmax + 1):(2nmax), (nmax + 1):(2nmax))
    𝐔 = zeros(CT, 2nmax, 2nmax)
    𝐔₁₁ = view(𝐔, 1:nmax, 1:nmax)
    𝐔₂₂ = view(𝐔, (nmax + 1):(2nmax), (nmax + 1):(2nmax))

    sc = Complex{R}(s)
    Mτd_q(n, n′) = q -> _sh_lookup(Mτd, n, n′, q, qlo, qhi)
    Mdτ_q(n, n′) = q -> _sh_lookup(Mdτ, n, n′, q, qlo, qhi)
    Mππττ_q(n) = q -> _sh_lookup(Mππττ, n, q, qlo, qhi)

    for n in 1:nmax, n′ in 1:nmax

        sym && isodd(n + n′) && continue
        if n != n′
            fτd = Mτd_q(n, n′)
            fdτ = Mdτ_q(n, n′)
            PL₁ = k * _sh_recon(ψ[n], ψ[n′], k, sc, fτd, 0)
            PL₂ = k * _sh_recon(ψ[n], ψ[n′], k, sc, fdτ, 0)
            PL₇ = k * _sh_recon(ψ′[n], ψ′[n′], k, sc, fτd, 0) +
                  (a[n] / sc / k) * _sh_recon(ψ[n], ψ[n′], k, sc, fτd, -2)
            PL₈ = k * _sh_recon(ψ′[n], ψ′[n′], k, sc, fdτ, 0) +
                  (a[n′] / sc / k) * _sh_recon(ψ[n], ψ[n′], k, sc, fdτ, -2)

            UL₁ = k * _sh_recon(χ[n], ψ[n′], k, sc, fτd, 0)
            UL₂ = k * _sh_recon(χ[n], ψ[n′], k, sc, fdτ, 0)
            UL₇ = k * _sh_recon(χ′[n], ψ′[n′], k, sc, fτd, 0) +
                  (a[n] / sc / k) * _sh_recon(χ[n], ψ[n′], k, sc, fτd, -2)
            UL₈ = k * _sh_recon(χ′[n], ψ′[n′], k, sc, fdτ, 0) +
                  (a[n′] / sc / k) * _sh_recon(χ[n], ψ[n′], k, sc, fdτ, -2)

            pref = 1im * A[n] * A[n′] * (sc^2 - 1) / (sc * (a[n] - a[n′]))
            𝐏₁₁[n, n′] = pref * (a[n] * PL₂ - a[n′] * PL₁)
            𝐏₂₂[n, n′] = pref * (a[n] * PL₈ - a[n′] * PL₇)
            𝐔₁₁[n, n′] = pref * (a[n] * UL₂ - a[n′] * UL₁)
            𝐔₂₂[n, n′] = pref * (a[n] * UL₈ - a[n′] * UL₇)
        else
            fππττ = Mππττ_q(n)
            fτd = Mτd_q(n, n)
            PL̃₁ = _sh_recon(ψ′[n], ψ[n], k, sc, fππττ, 0) -
                   sc * _sh_recon(ψ[n], ψ′[n], k, sc, fππττ, 0)
            PL̃₂ = sc * _sh_recon(ψ′[n], ψ[n], k, sc, fππττ, 0) -
                   _sh_recon(ψ[n], ψ′[n], k, sc, fππττ, 0)
            PL̃₃ = (1 / sc / k) * _sh_recon(ψ[n], ψ[n], k, sc, fτd, -2)

            UL̃₁ = _sh_recon(χ′[n], ψ[n], k, sc, fππττ, 0) -
                   sc * _sh_recon(χ[n], ψ′[n], k, sc, fππττ, 0)
            UL̃₂ = sc * _sh_recon(χ′[n], ψ[n], k, sc, fππττ, 0) -
                   _sh_recon(χ[n], ψ′[n], k, sc, fππττ, 0)
            UL̃₃ = (1 / sc / k) * _sh_recon(χ[n], ψ[n], k, sc, fτd, -2)

            𝐏₁₁[n, n] = -1im / sc * A[n]^2 * PL̃₁
            𝐏₂₂[n, n] = -1im / sc * A[n]^2 * (PL̃₂ + (sc^2 - 1) * a[n] * PL̃₃)
            𝐔₁₁[n, n] = -1im / sc * A[n]^2 * UL̃₁
            𝐔₂₂[n, n] = -1im / sc * A[n]^2 * (UL̃₂ + (sc^2 - 1) * a[n] * UL̃₃)
        end
    end

    return 𝐏, 𝐔
end
