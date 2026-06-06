# ── Two-layer solver API ──────────────────────────────────────────────────────
#
# The T-matrix of an axisymmetric scatterer can be obtained by several methods
# (classic EBCM, the cancellation-free `stable` EBCM, the Sh-matrix moment
# separation, IITM), each with its own discretization knobs. This file unifies
# them behind one verb, `transition_matrix(s, λ, solver)`, organized into two
# clearly separated layers:
#
#   • Fixed solution  — a solver carrying an explicit discretization; one build.
#       transition_matrix(s, λ, EBCM(nₘₐₓ, Ng; stable))
#       transition_matrix(s, λ, IITM(nₘₐₓ, Nr, Nϑ))
#       transition_matrix(s, λ, ShMatrix(nₘₐₓ, Ng))
#
#   • Iterative solver — wraps a method + convergence policy and sweeps the
#     discretization until the scattering/extinction efficiencies stabilize.
#       transition_matrix(s, λ, Iterative(EBCM; stable, threshold, …))
#
# Each fixed solver dispatches to the underlying implementation unchanged; the
# iterative driver reuses `routine_mishchenko` and the `m=0` efficiency probes.

"""
    AbstractSolver

Supertype of every T-matrix solver. A solver is either a fixed-discretization
solver ([`EBCM`](@ref), [`IITM`](@ref), [`ShMatrix`](@ref) — subtypes of
`AbstractFixedSolver`) or an [`Iterative`](@ref) wrapper that converges one.
"""
abstract type AbstractSolver end

"Supertype of fixed-discretization solvers (one build, no iteration)."
abstract type AbstractFixedSolver <: AbstractSolver end

# ── Fixed solvers ─────────────────────────────────────────────────────────────
"""
    EBCM(nₘₐₓ, Ng; stable = false)

Fixed-discretization Extended Boundary Condition Method solver: truncation order
`nₘₐₓ` and `Ng` Gauss–Legendre quadrature points. `stable = true` assembles the
`𝐔` matrix with the cancellation-free `F⁺` formulation (spheroids only). Used as
`transition_matrix(s, λ, EBCM(nₘₐₓ, Ng))`.
"""
struct EBCM <: AbstractFixedSolver
    nₘₐₓ::Int
    Ng::Int
    stable::Bool
end
EBCM(nₘₐₓ::Integer, Ng::Integer; stable::Bool = false) = EBCM(Int(nₘₐₓ), Int(Ng), stable)

"""
    IITM(nₘₐₓ, Nr, Nϑ; rₘᵢₙ = nothing)
    IITM(nₘₐₓ, Nr, Nϑ, Nφ; rₘᵢₙ = nothing)

Fixed-discretization invariant-imbedding (IITM) solver: order `nₘₐₓ`, `Nr` radial
and `Nϑ` zenithal quadrature points (plus `Nφ` azimuthal points for arbitrary /
N-fold shapes). `rₘᵢₙ` defaults to `rmin(s)` when left `nothing`.
"""
struct IITM{F} <: AbstractFixedSolver
    nₘₐₓ::Int
    Nr::Int
    Nϑ::Int
    Nφ::Union{Nothing, Int}
    rₘᵢₙ::F
end
function IITM(nₘₐₓ::Integer, Nr::Integer, Nϑ::Integer; rₘᵢₙ = nothing)
    IITM(Int(nₘₐₓ), Int(Nr), Int(Nϑ), nothing, rₘᵢₙ)
end
function IITM(nₘₐₓ::Integer, Nr::Integer, Nϑ::Integer, Nφ::Integer; rₘᵢₙ = nothing)
    IITM(Int(nₘₐₓ), Int(Nr), Int(Nϑ), Int(Nφ), rₘᵢₙ)
end

"""
    ShMatrix(nₘₐₓ, Ng; B = nothing, momtype = nothing, store = nothing)

Fixed-discretization Sh-matrix (moment-separation) solver. `B` is the number of
radial power-series terms; `momtype`/`store` control the moment precision (see
[`prepare_sh`](@ref); `nothing` selects the defaults). As a single build,
`transition_matrix(s, λ, ShMatrix(nₘₐₓ, Ng))`; for a parameter sweep,
`transition_matrix(s, λs, mᵣs, ShMatrix(nₘₐₓ, Ng))` prepares once and reuses.
"""
struct ShMatrix <: AbstractFixedSolver
    nₘₐₓ::Int
    Ng::Int
    B::Union{Nothing, Int}
    momtype::Union{Nothing, Type}
    store::Union{Nothing, Type}
end
function ShMatrix(nₘₐₓ::Integer, Ng::Integer; B = nothing, momtype = nothing,
        store = nothing)
    ShMatrix(Int(nₘₐₓ), Int(Ng), B, momtype, store)
end

# ── Fixed dispatch (thin shims over the existing implementations) ─────────────
function transition_matrix(s::AbstractAxisymmetricShape, λ, slv::EBCM)
    return transition_matrix(s, λ, slv.nₘₐₓ, slv.Ng; stable = slv.stable)
end

function transition_matrix(s::AbstractShape, λ, slv::IITM)
    rₘᵢₙ = slv.rₘᵢₙ === nothing ? rmin(s) : slv.rₘᵢₙ
    if slv.Nφ === nothing
        return transition_matrix_iitm(s, λ, slv.nₘₐₓ, slv.Nr, slv.Nϑ; rₘᵢₙ)
    else
        return transition_matrix_iitm(s, λ, slv.nₘₐₓ, slv.Nr, slv.Nϑ, slv.Nφ; rₘᵢₙ)
    end
end

function transition_matrix(s::AbstractAxisymmetricShape{T}, λ, slv::ShMatrix) where {T}
    B = slv.B === nothing ? max(30, slv.nₘₐₓ + 15) : slv.B
    store = slv.store === nothing ? T : slv.store
    prep = prepare_sh(s, slv.nₘₐₓ, slv.Ng; B, momtype = slv.momtype, store)
    return transition_matrix(prep, λ)
end

# Sh-matrix parameter sweep: prepare once, reconstruct at each (λ, mᵣ).
function transition_matrix(s::AbstractAxisymmetricShape{T}, λs::AbstractVector, mᵣs,
        slv::ShMatrix) where {T}
    B = slv.B === nothing ? max(30, slv.nₘₐₓ + 15) : slv.B
    store = slv.store === nothing ? T : slv.store
    prep = prepare_sh(s, slv.nₘₐₓ, slv.Ng; B, momtype = slv.momtype, store)
    return transition_matrix_spectrum(prep, λs, mᵣs)
end

# ── Iterative solver ──────────────────────────────────────────────────────────
@doc raw"""
    ConvergencePolicy(; threshold, ndgs, maxiter, nₘₐₓ_only, nₛₜₐᵣₜ, Ngₛₜₐᵣₜ, routine_generator)

Convergence settings for an [`Iterative`](@ref) solver. `threshold` is the
relative `Qsca`/`Qext` tolerance, `ndgs` the quadrature points added per order,
`nₛₜₐᵣₜ`/`Ngₛₜₐᵣₜ` the starting resolution (`0` ⇒ auto from ``k \cdot r_{\max}``),
`nₘₐₓ_only` stops once `nₘₐₓ` converges, and `routine_generator` builds the
stepping routine (default `routine_mishchenko`).
"""
struct ConvergencePolicy{G}
    threshold::Float64
    ndgs::Int
    maxiter::Int
    nₘₐₓ_only::Bool
    nₛₜₐᵣₜ::Int
    Ngₛₜₐᵣₜ::Int
    routine_generator::G
end
function ConvergencePolicy(; threshold = 1.0e-4, ndgs = 4, maxiter = 20,
        nₘₐₓ_only = false, nₛₜₐᵣₜ = 0, Ngₛₜₐᵣₜ = 0,
        routine_generator = routine_mishchenko)
    ConvergencePolicy(Float64(threshold), Int(ndgs), Int(maxiter), nₘₐₓ_only,
        Int(nₛₜₐᵣₜ), Int(Ngₛₜₐᵣₜ), routine_generator)
end

"""
    Iterative(EBCM; stable = false, threshold = 1e-4, ndgs = 4, maxiter = 20, …)

Iterative solver that sweeps the discretization of a fixed-solver method until
convergence. Currently the EBCM method is supported, including `stable = true`:

```julia
transition_matrix(s, λ, Iterative(EBCM))                  # classic, auto-converged
transition_matrix(s, λ, Iterative(EBCM; stable = true))   # stabilized + auto-converged
```

The non-resolution method options (e.g. `stable`) and the [`ConvergencePolicy`](@ref)
keywords are passed together; the method is named by its fixed-solver type.
"""
struct Iterative{F <: AbstractFixedSolver, O <: NamedTuple, P <: ConvergencePolicy} <:
       AbstractSolver
    opts::O
    policy::P
end

function Iterative(::Type{EBCM}; stable::Bool = false, kwargs...)
    opts = (; stable)
    policy = ConvergencePolicy(; kwargs...)
    return Iterative{EBCM, typeof(opts), typeof(policy)}(opts, policy)
end

# Generic driver: resolution-agnostic. Each method provides the four hooks below.
function transition_matrix(s::AbstractAxisymmetricShape, λ, it::Iterative)
    pol = it.policy
    routine = pol.routine_generator(pol.threshold, pol.ndgs, pol.nₘₐₓ_only)
    fixed = _solver_initial(it, s, λ)
    for _ in 1:(pol.maxiter)
        pr = _solver_probe(it, s, λ, fixed)
        nxt = _solver_step(it, fixed, pr.Qsca, pr.Qext, routine)
        nxt === nothing && return _solver_assemble(it, s, λ, fixed, pr.state)
        fixed = nxt
    end
    error("transition_matrix: failed to converge in $(pol.maxiter) iterations")
end

# ── Iterative interface — EBCM specialization (cheap m=0 probe) ───────────────
function _solver_initial(it::Iterative{EBCM}, s, λ)
    it.opts.stable && !(s isa Spheroid) &&
        throw(ArgumentError("Iterative(EBCM; stable=true) is only valid for spheroids: the \
              cancellation-free F⁺ integrands rely on the spheroid surface making \
              the divergent Laurent terms integrate to zero (Somerville et al. (2013))."))
    pol = it.policy
    nₛₜₐᵣₜ = pol.nₛₜₐᵣₜ
    if nₛₜₐᵣₜ == 0
        kr = 2π * rmax(s) / λ
        nₛₜₐᵣₜ = max(4, ceil(Int, kr + 4.05 * ∛kr))
    end
    Ng = pol.Ngₛₜₐᵣₜ == 0 ? nₛₜₐᵣₜ * pol.ndgs : pol.Ngₛₜₐᵣₜ
    return EBCM(nₛₜₐᵣₜ, Ng; stable = it.opts.stable)
end

function _solver_probe(it::Iterative{EBCM}, s, λ, f::EBCM)
    T₀,
    cache = transition_matrix_m₀(s, λ, f.nₘₐₓ, f.Ng; reuse = true,
        stable = it.opts.stable)
    return (; Qsca = scattering_efficiency_m₀(T₀), Qext = extinction_efficiency_m₀(T₀),
        state = (; T₀, cache))
end

function _solver_step(it::Iterative{EBCM}, f::EBCM, Qsca, Qext, routine)
    nₘₐₓ′, Ng′ = routine(f.nₘₐₓ, f.Ng, Qsca, Qext)
    return nₘₐₓ′ == -1 ? nothing : EBCM(nₘₐₓ′, Ng′; stable = f.stable)
end

function _solver_assemble(it::Iterative{EBCM}, s::AbstractAxisymmetricShape{T, CT}, λ,
        f::EBCM, state) where {T, CT}
    𝐓 = Vector{Matrix{CT}}(undef, f.nₘₐₓ + 1)
    𝐓[1] = state.T₀
    for m in 1:(f.nₘₐₓ)
        𝐓[m + 1] = transition_matrix_m(m, s, λ, f.nₘₐₓ, f.Ng; cache = state.cache,
            stable = it.opts.stable)
    end
    return AxisymmetricTransitionMatrix{CT, f.nₘₐₓ, typeof(𝐓), T}(𝐓)
end

# ── Iterative interface — generic fallback (extensibility seam) ───────────────
# Any iterable method without a cheap probe can converge via a full-T build and
# the decoupled cross-section probes. Not exercised in v1 (only EBCM is iterable);
# this is where `Iterative(IITM)` for axisymmetric shapes will plug in.
function _solver_probe(it::Iterative, s, λ, fixed::AbstractFixedSolver)
    𝐓 = transition_matrix(s, λ, fixed)
    return (; Qsca = scattering_cross_section(𝐓, λ), Qext = extinction_cross_section(𝐓, λ),
        state = (; 𝐓))
end
_solver_assemble(::Iterative, s, λ, ::AbstractFixedSolver, state) = state.𝐓

# ── Bare entry: no solver ⇒ auto-converged classic EBCM ───────────────────────
"""
    transition_matrix(s::AbstractAxisymmetricShape, λ)

Compute the T-matrix of an axisymmetric scatterer at wavelength `λ`, automatically
converging the truncation order and quadrature with classic EBCM. Equivalent to
`transition_matrix(s, λ, Iterative(EBCM))`. For a stabilized or fixed build, or
another method, pass an explicit solver (see [`AbstractSolver`](@ref)).
"""
function transition_matrix(s::AbstractAxisymmetricShape, λ)
    transition_matrix(s, λ, Iterative(EBCM))
end

@testitem "Fixed solver dispatch matches the underlying implementations" begin
    using TransitionMatrices: Spheroid, Cylinder, transition_matrix_iitm, prepare_sh,
                              calc_Csca
    λ = 2π
    s = Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)

    # EBCM fixed solver === the 4-arg call (same computation)
    @test transition_matrix(s, λ, EBCM(8, 64)).𝐓 == transition_matrix(s, λ, 8, 64).𝐓

    # EBCM stable fixed solver === the 4-arg stable call
    sp = Spheroid{Float64, ComplexF64}(2.5198421, 10.079368, 1.55 + 0.01im)
    @test transition_matrix(sp, λ, EBCM(20, 240; stable = true)).𝐓 ==
          transition_matrix(sp, λ, 20, 240; stable = true).𝐓

    # ShMatrix fixed solver === prepare_sh + assemble
    prep = prepare_sh(s, 8, 64)
    @test calc_Csca(transition_matrix(s, λ, ShMatrix(8, 64))) ≈
          calc_Csca(transition_matrix(prep, λ))

    # IITM fixed solver === transition_matrix_iitm
    c = Cylinder{Float64, ComplexF64}(1.0, 2.0, 1.5 + 0.02im)
    @test transition_matrix(c, λ, IITM(6, 16, 20)).𝐓 ==
          transition_matrix_iitm(c, λ, 6, 16, 20).𝐓
end

@testitem "Iterative(EBCM) auto-converges and equals the bare call" begin
    using TransitionMatrices: Spheroid, calc_Csca, calc_Cext
    λ = 2π
    s = Spheroid{Float64, ComplexF64}(2.0, 1.0, 1.5 + 0.02im)

    Ti = transition_matrix(s, λ, Iterative(EBCM))
    Tb = transition_matrix(s, λ)                       # bare === Iterative(EBCM)
    @test calc_Csca(Ti) == calc_Csca(Tb)

    Tf = transition_matrix(s, λ, EBCM(16, 200))        # high-resolution reference
    @test calc_Csca(Ti)≈calc_Csca(Tf) rtol=1e-3
    @test calc_Cext(Ti)≈calc_Cext(Tf) rtol=1e-3

    Tt = transition_matrix(s, λ, Iterative(EBCM; threshold = 1e-7))
    @test calc_Csca(Tt)≈calc_Csca(Tf) rtol=1e-5
end

@testitem "Iterative(EBCM; stable=true) converges a high-aspect spheroid" begin
    using TransitionMatrices: Spheroid, Cylinder, calc_Csca
    λ = 2π
    sp = Spheroid{Float64, ComplexF64}(2.5198421, 10.079368, 1.55 + 0.01im)

    Tst = transition_matrix(sp, λ, Iterative(EBCM; stable = true, threshold = 1e-6))
    Tref = transition_matrix(sp, λ, EBCM(28, 400; stable = true))
    @test calc_Csca(Tst)≈calc_Csca(Tref) rtol=1e-3

    # stable iterative is gated to spheroids
    c = Cylinder{Float64, ComplexF64}(1.0, 2.0, 1.5 + 0.02im)
    @test_throws ArgumentError transition_matrix(c, λ, Iterative(EBCM; stable = true))
end
