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

# ── Precomputed coefficient tables (size/material-independent) ────────────────
"""
    _sh_coeff_tables(nmax, B, R) -> NamedTuple(ψ, ψ′, χ, χ′)

Build the radial power-series coefficient tables once (they depend only on the
order and `B`, not on `k`/`mᵣ`), so a parameter sweep never re-runs the BigInt
`_β`/`_γ` arithmetic. Each entry is an `OffsetVector` over `n = 1:nmax` of
`(base, coeffs)` from [`_sh_series`](@ref).
"""
function _sh_coeff_tables(nmax::Integer, B::Integer, ::Type{R}) where {R}
    mk(kind) = OffsetArray([_sh_series(kind, n, B, R) for n in 1:nmax], 1:nmax)
    return (ψ = mk(:ψ), ψ′ = mk(:ψ′), χ = mk(:χ), χ′ = mk(:χ′))
end

# ── Per-evaluation power weighting ────────────────────────────────────────────
"""
    _sh_weighted(tbl, nmin, nmax, z) -> OffsetVector of (base, w)

Fold the argument power `z^{base+2b}` into the series coefficients, giving
`w[b+1] = coeffs[b+1]·z^{base+2b}`, with `z = k` for the regular factor `f(kr)`
or `z = s·k` for the internal factor `g(s·kr)`. Built once per `(k, mᵣ)` point
per order and reused across all integrands. Powers advance incrementally — no
`^` in the hot loop. The `base` is kept so [`_sh_conv`](@ref) can index the
shape moment by the correct `r`-power.
"""
function _sh_weighted(tbl, nmin::Integer, nmax::Integer, z::Number)
    z2 = z * z
    Tz = typeof(z * one(eltype(tbl[nmin][2])))
    out = OffsetArray(Vector{Tuple{Int, Vector{Tz}}}(undef, nmax - nmin + 1), nmin:nmax)
    for n in nmin:nmax
        base, cf = tbl[n]
        B = length(cf)
        w = Vector{Tz}(undef, B)
        zb = z^base
        @inbounds for b in 1:B
            w[b] = cf[b] * zb
            zb *= z2
        end
        out[n] = (base, w)
    end
    return out
end

# ── Folded coefficient×moment reconstruction ──────────────────────────────────
"""
    _sh_conv(f, g, M, shift) -> Complex

Reconstruct `∫ (geometry) · f_n(kr) · g_{n′}(s·kr) · r^{shift} dϑ` from the
weighted radial terms `f = (basef, u)`, `g = (baseg, v)` (see [`_sh_weighted`](@ref))
and the shape-moment lookup `M(q)`. The double sum is folded along the
anti-diagonals `j = b+c`: the `r`-power is `q₀ + 2j` (`q₀ = basef+baseg+shift`)
and depends only on `j`, so the moment `M` is fetched `O(B)` times instead of
`O(B²)`, and the inner convolution `Σ_{b+c=j} u[b]·v[c]` reuses the precomputed
weighted terms with no `^` and no per-term coefficient work.
"""
@inline function _sh_conv(f, g, M, shift::Integer)
    basef, u = f
    baseg, v = g
    Bf = length(u)
    Bg = length(v)
    q0 = basef + baseg + shift
    acc = zero(eltype(v))
    @inbounds for j in 0:(Bf + Bg - 2)
        si = zero(eltype(v))
        for b in max(0, j - (Bg - 1)):min(j, Bf - 1)
            si += u[b + 1] * v[j - b + 1]
        end
        acc += si * M(q0 + 2j)
    end
    return acc
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
    _sh_geometry(shape, Ng, qlo, qhi; momtype) -> NamedTuple

Compute the `m`-independent geometry once: the (high-precision) Gauss–Legendre
nodes/weights, the polar angles `ϑ`, and the two `r`-power weight tables shared
by every azimuthal block,

  WR[i,q] = wᵢ r′ᵢ rᵢ^q     (for the `τd`/`πd` families)
  Wq[i,q] = wᵢ rᵢ^q         (for the `ππττ` family, which carries no `r′`)

`gaussquad` and the `rᵢ^q` table do not depend on `m`, so hoisting them out of
the per-`m` loop avoids recomputing them `nmax+1` times.
"""
function _sh_geometry(shape, Ng::Integer, qlo::Integer, qhi::Integer;
        momtype::Type{R} = BigFloat) where {R}
    sym = has_symmetric_plane(shape)
    ng = sym ? Ng ÷ 2 : Ng

    sh = _sh_widen(shape, R)
    x, w, r, r′ = gaussquad(sh, Ng)
    ϑ = acos.(x)

    nq = qhi - qlo + 1
    WR = OffsetArray(zeros(R, ng, nq), 1:ng, qlo:qhi)
    Wq = OffsetArray(zeros(R, ng, nq), 1:ng, qlo:qhi)
    for i in 1:ng
        ri = r[i]
        wi = w[i]
        wri = wi * r′[i]
        rqᵢ = one(R)
        Wq[i, 0] = wi
        WR[i, 0] = wri
        for q in 1:qhi
            rqᵢ *= ri
            Wq[i, q] = wi * rqᵢ
            WR[i, q] = wri * rqᵢ
        end
        rqᵢ = one(R)
        for q in -1:-1:qlo
            rqᵢ /= ri
            Wq[i, q] = wi * rqᵢ
            WR[i, q] = wri * rqᵢ
        end
    end

    return (; sym, ng, ϑ, WR, Wq, qlo, qhi)
end

"""
    _sh_moments_m(geom, m, nmax; store) -> NamedTuple

Compute the shape-only moment tables for azimuthal index `m` from the shared
[`_sh_geometry`](@ref) `geom`, accumulated at the geometry's precision and stored
as `store` (default `Float64`). Families (sum over `i in 1:ng`, mirroring
`ebcm_matrices_m₀`):

  Mτd[n,n′,q]  = Σ WR[i,q] τ[i,n] d[i,n′]
  Mπd[n,n′,q]  = Σ WR[i,q] π[i,n] d[i,n′]   (only `m>0`)
  Mππττ[n,q]   = Σ Wq[i,q] (π[i,n]²+τ[i,n]²)

`Mdτ[n,n′,q] = Mτd[n′,n,q]` exactly (swap the two orders), so it is not stored —
the assembly reads `Mτd` with the indices transposed. The angular products are
formed once per `(n,n′)` instead of once per `q`.
"""
function _sh_moments_m(geom, m::Integer, nmax::Integer;
        store::Type{Ts} = Float64) where {Ts}
    WR, Wq = geom.WR, geom.Wq
    R = eltype(WR)
    ϑ = geom.ϑ
    ng = geom.ng
    qlo, qhi = geom.qlo, geom.qhi
    Ng = length(ϑ)
    mm = Int(m)
    nmin = max(1, mm)
    nq = qhi - qlo + 1

    d = OffsetArray(zeros(R, Ng, nmax - mm + 1), 1:Ng, mm:nmax)
    τ = similar(d)
    𝜋 = similar(d)
    for i in 1:Ng
        wigner_d_recursion!(view(d, i, :), 0, mm, nmax, ϑ[i]; deriv = view(τ, i, :))
        for n in mm:nmax
            𝜋[i, n] = pi_func(R, mm, n, ϑ[i]; d = d[i, n])
        end
    end

    Mτd = OffsetArray(zeros(Ts, nmax, nmax, nq), 1:nmax, 1:nmax, qlo:qhi)
    Mπd = OffsetArray(zeros(Ts, nmax, nmax, nq), 1:nmax, 1:nmax, qlo:qhi)
    Mππττ = OffsetArray(zeros(Ts, nmax, nq), 1:nmax, qlo:qhi)

    gτd = Vector{R}(undef, ng)
    gπd = Vector{R}(undef, ng)
    for n in nmin:nmax, n′ in nmin:nmax

        @inbounds for i in 1:ng
            dᵢ = d[i, n′]
            gτd[i] = τ[i, n] * dᵢ
            gπd[i] = 𝜋[i, n] * dᵢ
        end
        @inbounds for q in qlo:qhi
            sτd = zero(R)
            sπd = zero(R)
            for i in 1:ng
                wr = WR[i, q]
                sτd += wr * gτd[i]
                sπd += wr * gπd[i]
            end
            Mτd[n, n′, q] = Ts(sτd)
            Mπd[n, n′, q] = Ts(sπd)
        end
    end

    gππττ = Vector{R}(undef, ng)
    for n in nmin:nmax
        @inbounds for i in 1:ng
            gππττ[i] = 𝜋[i, n]^2 + τ[i, n]^2
        end
        @inbounds for q in qlo:qhi
            sm = zero(R)
            for i in 1:ng
                sm += Wq[i, q] * gππττ[i]
            end
            Mππττ[n, q] = Ts(sm)
        end
    end

    return (; Mτd, Mπd, Mππττ, nmin, qlo, qhi, sym = geom.sym, ng)
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
    _sh_matrices_m₀(mom, coeffs, k, s, nmax, CT) -> (𝐏, 𝐔)

Reconstruct the `m=0` `𝐏` and `𝐔` blocks from the precomputed `m=0` moment
tables `mom` and the size/material-independent coefficient tables `coeffs` at
size parameter `k` and refractive index `s`. The algebra reproduces
`ebcm_matrices_m₀` term by term.
"""
function _sh_matrices_m₀(mom, coeffs, k::Real, s::Number, nmax::Integer,
        ::Type{CT}) where {CT}
    R = real(CT)
    sym = mom.sym
    qlo, qhi = mom.qlo, mom.qhi
    Mτd, Mππττ = mom.Mτd, mom.Mππττ

    a = [n * (n + 1) for n in 1:nmax]
    A = [√(R(2n + 1) / (2n * (n + 1))) for n in 1:nmax]

    sc = Complex{R}(s)
    kk = R(k)
    sk = sc * kk
    Uψ = _sh_weighted(coeffs.ψ, 1, nmax, kk)
    Uψ′ = _sh_weighted(coeffs.ψ′, 1, nmax, kk)
    Uχ = _sh_weighted(coeffs.χ, 1, nmax, kk)
    Uχ′ = _sh_weighted(coeffs.χ′, 1, nmax, kk)
    Vψ = _sh_weighted(coeffs.ψ, 1, nmax, sk)
    Vψ′ = _sh_weighted(coeffs.ψ′, 1, nmax, sk)

    𝐏 = zeros(CT, 2nmax, 2nmax)
    𝐏₁₁ = view(𝐏, 1:nmax, 1:nmax)
    𝐏₂₂ = view(𝐏, (nmax + 1):(2nmax), (nmax + 1):(2nmax))
    𝐔 = zeros(CT, 2nmax, 2nmax)
    𝐔₁₁ = view(𝐔, 1:nmax, 1:nmax)
    𝐔₂₂ = view(𝐔, (nmax + 1):(2nmax), (nmax + 1):(2nmax))

    # Mdτ[n,n′,q] = Mτd[n′,n,q] exactly, so read Mτd with the orders swapped.
    Mτd_q(n, n′) = q -> _sh_lookup(Mτd, n, n′, q, qlo, qhi)
    Mdτ_q(n, n′) = q -> _sh_lookup(Mτd, n′, n, q, qlo, qhi)
    Mππττ_q(n) = q -> _sh_lookup(Mππττ, n, q, qlo, qhi)

    for n in 1:nmax, n′ in 1:nmax

        sym && isodd(n + n′) && continue
        if n != n′
            fτd = Mτd_q(n, n′)
            fdτ = Mdτ_q(n, n′)
            PL₁ = kk * _sh_conv(Uψ[n], Vψ[n′], fτd, 0)
            PL₂ = kk * _sh_conv(Uψ[n], Vψ[n′], fdτ, 0)
            PL₇ = kk * _sh_conv(Uψ′[n], Vψ′[n′], fτd, 0) +
                  (a[n] / sc / kk) * _sh_conv(Uψ[n], Vψ[n′], fτd, -2)
            PL₈ = kk * _sh_conv(Uψ′[n], Vψ′[n′], fdτ, 0) +
                  (a[n′] / sc / kk) * _sh_conv(Uψ[n], Vψ[n′], fdτ, -2)

            UL₁ = kk * _sh_conv(Uχ[n], Vψ[n′], fτd, 0)
            UL₂ = kk * _sh_conv(Uχ[n], Vψ[n′], fdτ, 0)
            UL₇ = kk * _sh_conv(Uχ′[n], Vψ′[n′], fτd, 0) +
                  (a[n] / sc / kk) * _sh_conv(Uχ[n], Vψ[n′], fτd, -2)
            UL₈ = kk * _sh_conv(Uχ′[n], Vψ′[n′], fdτ, 0) +
                  (a[n′] / sc / kk) * _sh_conv(Uχ[n], Vψ[n′], fdτ, -2)

            pref = 1im * A[n] * A[n′] * (sc^2 - 1) / (sc * (a[n] - a[n′]))
            𝐏₁₁[n, n′] = pref * (a[n] * PL₂ - a[n′] * PL₁)
            𝐏₂₂[n, n′] = pref * (a[n] * PL₈ - a[n′] * PL₇)
            𝐔₁₁[n, n′] = pref * (a[n] * UL₂ - a[n′] * UL₁)
            𝐔₂₂[n, n′] = pref * (a[n] * UL₈ - a[n′] * UL₇)
        else
            fππττ = Mππττ_q(n)
            fτd = Mτd_q(n, n)
            PL̃₁ = _sh_conv(Uψ′[n], Vψ[n], fππττ, 0) -
                   sc * _sh_conv(Uψ[n], Vψ′[n], fππττ, 0)
            PL̃₂ = sc * _sh_conv(Uψ′[n], Vψ[n], fππττ, 0) -
                   _sh_conv(Uψ[n], Vψ′[n], fππττ, 0)
            PL̃₃ = (1 / sc / kk) * _sh_conv(Uψ[n], Vψ[n], fτd, -2)

            UL̃₁ = _sh_conv(Uχ′[n], Vψ[n], fππττ, 0) -
                   sc * _sh_conv(Uχ[n], Vψ′[n], fππττ, 0)
            UL̃₂ = sc * _sh_conv(Uχ′[n], Vψ[n], fππττ, 0) -
                   _sh_conv(Uχ[n], Vψ′[n], fππττ, 0)
            UL̃₃ = (1 / sc / kk) * _sh_conv(Uχ[n], Vψ[n], fτd, -2)

            𝐏₁₁[n, n] = -1im / sc * A[n]^2 * PL̃₁
            𝐏₂₂[n, n] = -1im / sc * A[n]^2 * (PL̃₂ + (sc^2 - 1) * a[n] * PL̃₃)
            𝐔₁₁[n, n] = -1im / sc * A[n]^2 * UL̃₁
            𝐔₂₂[n, n] = -1im / sc * A[n]^2 * (UL̃₂ + (sc^2 - 1) * a[n] * UL̃₃)
        end
    end

    return 𝐏, 𝐔
end

# ── m>0 assembly from moments (mirrors `_transition_matrix_m_core`) ───────────
"""
    _sh_matrices_m(mom, coeffs, m, k, s, nmax, CT) -> (𝐏, 𝐔)

Reconstruct the `m`-th `𝐏` and `𝐔` blocks (`2nn × 2nn`, `nn = nmax-m+1` for
`m≥1`) from the precomputed moment tables `mom` for that `m` and the
size/material-independent coefficient tables `coeffs`. The `K`-blocks use the
`M^{πd}` family; the `L`-blocks reuse the same reconstruction as the `m=0` case
(with the `m`-dependent Wigner functions baked into `mom`). The algebra
reproduces `_transition_matrix_m_core` term by term.
"""
function _sh_matrices_m(mom, coeffs, m::Integer, k::Real, s::Number, nmax::Integer,
        ::Type{CT}) where {CT}
    R = real(CT)
    sym = mom.sym
    qlo, qhi = mom.qlo, mom.qhi
    Mτd, Mπd, Mππττ = mom.Mτd, mom.Mπd, mom.Mππττ
    nₘᵢₙ = max(1, Int(m))
    nn = nmax - nₘᵢₙ + 1

    a = OffsetArray([n * (n + 1) for n in nₘᵢₙ:nmax], nₘᵢₙ:nmax)
    A = OffsetArray([√(R(2n + 1) / (2n * (n + 1))) for n in nₘᵢₙ:nmax], nₘᵢₙ:nmax)

    sc = Complex{R}(s)
    kk = R(k)
    sk = sc * kk
    Uψ = _sh_weighted(coeffs.ψ, nₘᵢₙ, nmax, kk)
    Uψ′ = _sh_weighted(coeffs.ψ′, nₘᵢₙ, nmax, kk)
    Uχ = _sh_weighted(coeffs.χ, nₘᵢₙ, nmax, kk)
    Uχ′ = _sh_weighted(coeffs.χ′, nₘᵢₙ, nmax, kk)
    Vψ = _sh_weighted(coeffs.ψ, nₘᵢₙ, nmax, sk)
    Vψ′ = _sh_weighted(coeffs.ψ′, nₘᵢₙ, nmax, sk)

    𝐏 = zeros(CT, 2nn, 2nn)
    𝐏₁₁ = OffsetArray(view(𝐏, 1:nn, 1:nn), nₘᵢₙ:nmax, nₘᵢₙ:nmax)
    𝐏₁₂ = OffsetArray(view(𝐏, 1:nn, (nn + 1):(2nn)), nₘᵢₙ:nmax, nₘᵢₙ:nmax)
    𝐏₂₁ = OffsetArray(view(𝐏, (nn + 1):(2nn), 1:nn), nₘᵢₙ:nmax, nₘᵢₙ:nmax)
    𝐏₂₂ = OffsetArray(view(𝐏, (nn + 1):(2nn), (nn + 1):(2nn)), nₘᵢₙ:nmax, nₘᵢₙ:nmax)

    𝐔 = zeros(CT, 2nn, 2nn)
    𝐔₁₁ = OffsetArray(view(𝐔, 1:nn, 1:nn), nₘᵢₙ:nmax, nₘᵢₙ:nmax)
    𝐔₁₂ = OffsetArray(view(𝐔, 1:nn, (nn + 1):(2nn)), nₘᵢₙ:nmax, nₘᵢₙ:nmax)
    𝐔₂₁ = OffsetArray(view(𝐔, (nn + 1):(2nn), 1:nn), nₘᵢₙ:nmax, nₘᵢₙ:nmax)
    𝐔₂₂ = OffsetArray(view(𝐔, (nn + 1):(2nn), (nn + 1):(2nn)), nₘᵢₙ:nmax, nₘᵢₙ:nmax)

    # Mdτ[n,n′,q] = Mτd[n′,n,q] exactly, so read Mτd with the orders swapped.
    Mτd_q(n, n′) = q -> _sh_lookup(Mτd, n, n′, q, qlo, qhi)
    Mdτ_q(n, n′) = q -> _sh_lookup(Mτd, n′, n, q, qlo, qhi)
    Mπd_q(n, n′) = q -> _sh_lookup(Mπd, n, n′, q, qlo, qhi)
    Mππττ_q(n) = q -> _sh_lookup(Mππττ, n, q, qlo, qhi)

    for n in nₘᵢₙ:nmax, n′ in nₘᵢₙ:nmax

        if !(sym && iseven(n + n′))
            fπd = Mπd_q(n, n′)
            PK₁ = kk * _sh_conv(Uψ[n], Vψ′[n′], fπd, 0)
            PK₂ = kk * _sh_conv(Uψ′[n], Vψ[n′], fπd, 0)
            UK₁ = kk * _sh_conv(Uχ[n], Vψ′[n′], fπd, 0)
            UK₂ = kk * _sh_conv(Uχ′[n], Vψ[n′], fπd, 0)

            𝐏₁₂[n, n′] = A[n] * A[n′] * (sc^2 - 1) / sc * PK₁
            𝐏₂₁[n, n′] = A[n] * A[n′] * (1 - sc^2) / sc * PK₂
            𝐔₁₂[n, n′] = A[n] * A[n′] * (sc^2 - 1) / sc * UK₁
            𝐔₂₁[n, n′] = A[n] * A[n′] * (1 - sc^2) / sc * UK₂
        end

        if !(sym && isodd(n + n′))
            if n != n′
                fτd = Mτd_q(n, n′)
                fdτ = Mdτ_q(n, n′)
                PL₁ = kk * _sh_conv(Uψ[n], Vψ[n′], fτd, 0)
                PL₂ = kk * _sh_conv(Uψ[n], Vψ[n′], fdτ, 0)
                PL₇ = kk * _sh_conv(Uψ′[n], Vψ′[n′], fτd, 0) +
                      (a[n] / sc / kk) * _sh_conv(Uψ[n], Vψ[n′], fτd, -2)
                PL₈ = kk * _sh_conv(Uψ′[n], Vψ′[n′], fdτ, 0) +
                      (a[n′] / sc / kk) * _sh_conv(Uψ[n], Vψ[n′], fdτ, -2)

                UL₁ = kk * _sh_conv(Uχ[n], Vψ[n′], fτd, 0)
                UL₂ = kk * _sh_conv(Uχ[n], Vψ[n′], fdτ, 0)
                UL₇ = kk * _sh_conv(Uχ′[n], Vψ′[n′], fτd, 0) +
                      (a[n] / sc / kk) * _sh_conv(Uχ[n], Vψ[n′], fτd, -2)
                UL₈ = kk * _sh_conv(Uχ′[n], Vψ′[n′], fdτ, 0) +
                      (a[n′] / sc / kk) * _sh_conv(Uχ[n], Vψ[n′], fdτ, -2)

                pref = 1im * A[n] * A[n′] * (sc^2 - 1) / (sc * (a[n] - a[n′]))
                𝐏₁₁[n, n′] = pref * (a[n] * PL₂ - a[n′] * PL₁)
                𝐏₂₂[n, n′] = pref * (a[n] * PL₈ - a[n′] * PL₇)
                𝐔₁₁[n, n′] = pref * (a[n] * UL₂ - a[n′] * UL₁)
                𝐔₂₂[n, n′] = pref * (a[n] * UL₈ - a[n′] * UL₇)
            else
                fππττ = Mππττ_q(n)
                fτd = Mτd_q(n, n)
                PL̃₁ = _sh_conv(Uψ′[n], Vψ[n], fππττ, 0) -
                       sc * _sh_conv(Uψ[n], Vψ′[n], fππττ, 0)
                PL̃₂ = sc * _sh_conv(Uψ′[n], Vψ[n], fππττ, 0) -
                       _sh_conv(Uψ[n], Vψ′[n], fππττ, 0)
                PL̃₃ = (1 / sc / kk) * _sh_conv(Uψ[n], Vψ[n], fτd, -2)

                UL̃₁ = _sh_conv(Uχ′[n], Vψ[n], fππττ, 0) -
                       sc * _sh_conv(Uχ[n], Vψ′[n], fππττ, 0)
                UL̃₂ = sc * _sh_conv(Uχ′[n], Vψ[n], fππττ, 0) -
                       _sh_conv(Uχ[n], Vψ′[n], fππττ, 0)
                UL̃₃ = (1 / sc / kk) * _sh_conv(Uχ[n], Vψ[n], fτd, -2)

                𝐏₁₁[n, n] = -1im / sc * A[n]^2 * PL̃₁
                𝐏₂₂[n, n] = -1im / sc * A[n]^2 * (PL̃₂ + (sc^2 - 1) * a[n] * PL̃₃)
                𝐔₁₁[n, n] = -1im / sc * A[n]^2 * UL̃₁
                𝐔₂₂[n, n] = -1im / sc * A[n]^2 * (UL̃₂ + (sc^2 - 1) * a[n] * UL̃₃)
            end
        end
    end

    return 𝐏, 𝐔
end

# ── Two-step public API: `prepare_sh` then `transition_matrix(prep, λ, mᵣ)` ───
"""
    ShPreparation

Precomputed shape-only moment tables for the Sh-matrix moment-separation method,
returned by [`prepare_sh`](@ref). Holds the moments for every azimuthal index
`m = 0:nmax`, which depend only on the particle geometry — not on the wavelength
or refractive index — so a single preparation feeds an arbitrarily long sweep of
`transition_matrix(prep, λ, mᵣ)` evaluations cheaply.

The `stable` field records whether the analytic `𝐔` stabilization is valid: it
is `true` only for spheroids, whose negative-`r`-power moments vanish (computed
at high precision they contribute ≈0, removing the irregular-product
cancellation). For other axisymmetric shapes the moment machinery and `𝐏` are
still correct, but the `𝐔` reconstruction is not stabilized.
"""
struct ShPreparation{ST, MT, CF}
    shape::ST
    nmax::Int
    Ng::Int
    B::Int
    moments::MT
    coeffs::CF
    stable::Bool
end

"""
    prepare_sh(shape, nmax, Ng; B, momtype, store) -> ShPreparation

Precompute the Sh-matrix shape-only moment tables for `shape` up to order `nmax`
with `Ng` Gauss–Legendre points. The result is reused across many
`transition_matrix(prep, λ, mᵣ)` calls (a wavelength / refractive-index sweep)
at the cost of only a cheap coefficient×moment sum each.

Keyword arguments:

- `B`: number of radial power-series terms per Riccati–Bessel function. It sets
  the largest size parameter the reconstruction resolves (the series, like the
  `F⁺` evaluation, needs more terms as `k·rₘₐₓ·|mᵣ|` grows); it is fixed at
  preparation time because it determines the moment band. Default `max(30,
  nmax+15)`.
- `momtype`: precision used to accumulate the moments (default `BigFloat`, so a
  spheroid's vanishing negative-power moments are accurate enough to cancel the
  irregular-product blow-up).
- `store`: element type the moments are stored as (default the shape's real
  type); reconstruction runs in this precision.

The analytic `𝐔` stabilization is enabled only for `Spheroid` (see
[`ShPreparation`](@ref)).
"""
function prepare_sh(shape::AbstractAxisymmetricShape{T}, nmax::Integer, Ng::Integer;
        B::Integer = max(30, nmax + 15), momtype::Type = BigFloat,
        store::Type = T) where {T}
    @assert iseven(Ng) "Ng must be even!"
    qlo, qhi = _sh_qband(nmax, B)
    geom = _sh_geometry(shape, Ng, qlo, qhi; momtype)
    moments = [_sh_moments_m(geom, m, nmax; store) for m in 0:nmax]
    coeffs = _sh_coeff_tables(nmax, B, store)
    return ShPreparation(shape, Int(nmax), Int(Ng), Int(B), moments, coeffs,
        shape isa Spheroid)
end

"""
    transition_matrix(prep::ShPreparation, λ, mᵣ) -> AxisymmetricTransitionMatrix
    transition_matrix(prep::ShPreparation, λ)      # uses the prepared shape's mᵣ

Assemble the full T-matrix at wavelength `λ` and relative refractive index `mᵣ`
from a [`prepare_sh`](@ref) preparation, reconstructing every `m`-block from the
precomputed moments. This is the cheap inner call of a parameter sweep; the
expensive geometry quadrature was done once in `prepare_sh`.
"""
function transition_matrix(prep::ShPreparation, λ::Real, mᵣ::Number)
    nmax = prep.nmax
    R = promote_type(eltype(prep.moments[1].Mτd), typeof(float(λ)),
        typeof(float(real(mᵣ))))
    CT = Complex{R}
    k = 2R(π) / R(λ)
    s = CT(mᵣ)
    coeffs = prep.coeffs

    𝐓 = Vector{Matrix{CT}}(undef, nmax + 1)
    𝐏, 𝐔 = _sh_matrices_m₀(prep.moments[1], coeffs, k, s, nmax, CT)
    𝐓[1] = 𝐓_from_𝐏_and_𝐔(𝐏, 𝐔)
    for m in 1:nmax
        𝐏ₘ, 𝐔ₘ = _sh_matrices_m(prep.moments[m + 1], coeffs, m, k, s, nmax, CT)
        𝐓[m + 1] = 𝐓_from_𝐏_and_𝐔(𝐏ₘ, 𝐔ₘ)
    end

    return AxisymmetricTransitionMatrix{CT, nmax, typeof(𝐓), R}(𝐓)
end

transition_matrix(prep::ShPreparation, λ::Real) = transition_matrix(prep, λ, prep.shape.m)

"""
    transition_matrix_spectrum(prep::ShPreparation, λs, mᵣs) -> Vector

Reconstruct a T-matrix at each `(λ, mᵣ)` pair, reusing the prepared moments. If
`mᵣs` is a single number it is held fixed across all `λs`; otherwise `λs` and
`mᵣs` are iterated together (and must have equal length). Pass a dispersion
table as `mᵣs` to sweep a material with wavelength-dependent index.
"""
function transition_matrix_spectrum(prep::ShPreparation, λs, mᵣs)
    ms = mᵣs isa Number ? Iterators.repeated(mᵣs, length(λs)) : mᵣs
    return [transition_matrix(prep, λ, m) for (λ, m) in zip(λs, ms)]
end

@testitem "Sh-matrix reproduces the ebcm P and U blocks (m=0 and m>0)" begin
    using TransitionMatrices: Spheroid, ebcm_matrices_m₀, ebcm_matrices_m,
                              _sh_geometry, _sh_moments_m, _sh_matrices_m₀,
                              _sh_matrices_m, _sh_qband, _sh_coeff_tables

    relmax(A,
        B;
        tol = 1e-10) = maximum(
        (abs(B[i, j]) < tol ? 0.0 : abs(A[i, j] - B[i, j]) / abs(B[i, j])
        for i in axes(A, 1), j in axes(A, 2)); init = 0.0)

    a, c, λ, sm = 2.0, 1.0, 2π, 1.5 + 0.02im
    nmax, Ng, B = 6, 120, 26
    shape = Spheroid{Float64, ComplexF64}(a, c, sm)
    k = 2π / λ
    qlo, qhi = _sh_qband(nmax, B)
    coeffs = _sh_coeff_tables(nmax, B, Float64)
    geom = _sh_geometry(shape, Ng, qlo, qhi)

    mom0 = _sh_moments_m(geom, 0, nmax)
    Psh, Ush = _sh_matrices_m₀(mom0, coeffs, k, sm, nmax, ComplexF64)
    Pref, Uref = ebcm_matrices_m₀(shape, λ, nmax, Ng)
    @test relmax(Psh, Pref) < 1e-9
    @test relmax(Ush, Uref) < 1e-7

    for m in 1:4
        mom = _sh_moments_m(geom, m, nmax)
        Pm, Um = _sh_matrices_m(mom, coeffs, m, k, sm, nmax, ComplexF64)
        Pr, Ur = ebcm_matrices_m(m, shape, λ, nmax, Ng)
        @test relmax(Pm, Pr) < 1e-9
        @test relmax(Um, Ur) < 1e-7
    end
end

@testitem "Sh-matrix transition_matrix reproduces the package T and cross sections" begin
    using TransitionMatrices: Spheroid, prepare_sh, transition_matrix, calc_Csca, calc_Cext

    # Moderate aspect-2 spheroid: full T must match the standard assembly.
    sh = Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)
    nmax, Ng, λ = 8, 160, 2π
    Tpkg = transition_matrix(sh, λ, nmax, Ng)
    prep = prepare_sh(sh, nmax, Ng; B = 28)
    Tsh = transition_matrix(prep, λ)
    @test prep.stable
    @test calc_Csca(Tsh) ≈ calc_Csca(Tpkg) rtol = 1e-9
    @test calc_Cext(Tsh) ≈ calc_Cext(Tpkg) rtol = 1e-9

    # High-aspect prolate: matches the stabilized (F⁺) path where the standard
    # Float64 assembly loses precision.
    shp = Spheroid{Float64, ComplexF64}(2.5198421, 10.079368, 1.55 + 0.01im)
    nmax2, Ng2 = 22, 360
    Tstab = transition_matrix(shp, λ, nmax2, Ng2; stable = true)
    prep2 = prepare_sh(shp, nmax2, Ng2; B = 42)
    Tsh2 = transition_matrix(prep2, λ)
    @test calc_Csca(Tsh2) ≈ calc_Csca(Tstab) rtol = 1e-6
end

@testitem "Sh-matrix spectrum sweep matches per-point assembly" begin
    using TransitionMatrices: Spheroid, prepare_sh, transition_matrix,
                              transition_matrix_spectrum, calc_Csca

    sh = Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)
    nmax, Ng = 8, 160
    prep = prepare_sh(sh, nmax, Ng; B = 28)

    λs = [2π, 2.3π, 2.6π]
    ms = [1.5 + 0.02im, 1.52 + 0.02im, 1.55 + 0.03im]
    Ts = transition_matrix_spectrum(prep, λs, ms)
    @test length(Ts) == length(λs)
    for (Ti, λi, mi) in zip(Ts, λs, ms)
        Tp = transition_matrix(Spheroid{Float64, ComplexF64}(2.0, 1.0, mi), λi, nmax, Ng)
        @test calc_Csca(Ti) ≈ calc_Csca(Tp) rtol = 1e-9
    end

    # A scalar refractive index is held fixed across the whole sweep.
    Tsf = transition_matrix_spectrum(prep, λs, sh.m)
    @test length(Tsf) == length(λs)
    @test calc_Csca(Tsf[1]) ≈ calc_Csca(transition_matrix(prep, λs[1], sh.m)) rtol = 1e-12
end
