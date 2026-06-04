# Benchmark suite for TransitionMatrices.jl.
#
# This file follows the de-facto Julia convention: it defines a single global
# `const SUITE::BenchmarkGroup`. The same file is consumed by both
# PkgBenchmark.jl (`benchmarkpkg`, `judge`) and AirspeedVelocity.jl
# (`benchpkg`), so no tooling lock-in is introduced.
#
# Run locally with `just bench`, or directly:
#   julia --project=benchmark -e 'include("benchmark/runbenchmarks.jl")'
#
# Groups
#   special_functions  – Riccati-Bessel and Wigner-d recursions (inner kernels)
#   ebcm               – EBCM blocks, fixed-order, and auto-converging solves
#   iitm               – Invariant Imbedding T-Matrix (axisymmetric + N-fold)
#   postprocessing     – far-field observables from a precomputed T-matrix
#   linearization      – analytical Jacobian (baseline for the inv→lu work)
#   precision          – Float64 vs Double64 on the same EBCM block

using BenchmarkTools
using TransitionMatrices

const SUITE = BenchmarkGroup()

# --------------------------------------------------------------------------
# Representative scatterers. With λ = 2π the size parameter is x = k·rmax = rmax.
# A mildly absorbing refractive index (1.5 + 0.02im) keeps the problems
# well-conditioned in Float64 while exercising the complex paths.
# --------------------------------------------------------------------------
const λ = 2π
const SPHEROID = Spheroid(2.0, 1.0, complex(1.5, 0.02))
const CYLINDER = Cylinder(1.0, 2.0, complex(1.5, 0.02))
const CHEBYSHEV = Chebyshev(1.0, 0.1, 4, complex(1.5, 0.02))
const PRISM = Prism(6, 1.0, 2.0, complex(1.5, 0.02))

# A reference T-matrix shared by the post-processing benchmarks. Built once at
# suite-load time so its cost stays out of the per-sample timings.
const T_REF = transition_matrix(SPHEROID, λ, 8, 32)
const θs = collect(range(0, π; length = 181))

# --------------------------------------------------------------------------
# Special functions — the innermost recursions, sized by truncation order.
# --------------------------------------------------------------------------
let g = SUITE["special_functions"] = BenchmarkGroup(["math", "recurrence"])
    for nmax in (20, 50, 100)
        x = float(nmax)  # argument near the order: the demanding regime
        nextra = TransitionMatrices.estimate_ricattibesselj_extra_terms(nmax, x)
        g["ricattibesselj_n$nmax"] = @benchmarkable TransitionMatrices.ricattibesselj($nmax,
                                                                                      $nextra,
                                                                                      $x)
    end
    for smax in (20, 50, 100)
        g["wigner_d_recursion_s$smax"] = @benchmarkable TransitionMatrices.wigner_d_recursion(0,
                                                                                              2,
                                                                                              $smax,
                                                                                              0.7;
                                                                                              deriv = true)
    end
end

# --------------------------------------------------------------------------
# EBCM — block construction, fixed-order solves, and end-to-end convergence.
# --------------------------------------------------------------------------
let g = SUITE["ebcm"] = BenchmarkGroup(["core", "linalg"])
    # The m = 0 block is the per-iteration hot path of the convergence loop.
    for nmax in (8, 12, 16)
        Ng = 4nmax
        g["m0_block_n$nmax"] = @benchmarkable transition_matrix_m₀($SPHEROID, $λ, $nmax, $Ng)
    end

    # Full fixed-order T-matrix across the three axisymmetric shape families.
    g["full_spheroid"] = @benchmarkable transition_matrix($SPHEROID, $λ, 8, 32)
    g["full_cylinder"] = @benchmarkable transition_matrix($CYLINDER, $λ, 8, 32)
    g["full_chebyshev"] = @benchmarkable transition_matrix($CHEBYSHEV, $λ, 8, 32)

    # Realistic user call: automatic nmax/Ng convergence from scratch.
    g["auto_converge_spheroid"] = @benchmarkable calc_T($SPHEROID, $λ) seconds=30
end

# --------------------------------------------------------------------------
# IITM — the numerically stable solver; the GPU/N-fold work targets this path.
# --------------------------------------------------------------------------
let g = SUITE["iitm"] = BenchmarkGroup(["core", "linalg"])
    g["axisym_spheroid"] = @benchmarkable calc_T_iitm($SPHEROID, $λ, 6, 20, 24) seconds=60
    g["nfold_prism6"] = @benchmarkable calc_T_iitm($PRISM, $λ, 6, 10, 24, 24) seconds=60
end

# --------------------------------------------------------------------------
# Post-processing — far-field observables from a precomputed T-matrix.
# The analytic vs numerical orientation average contrast is intentional: it
# documents the (already implemented) speedup of Mishchenko's closed form.
# --------------------------------------------------------------------------
let g = SUITE["postprocessing"] = BenchmarkGroup(["farfield"])
    g["amplitude_matrix"] = @benchmarkable amplitude_matrix($T_REF, 0.0, 0.0, π / 4, 0.0;
                                                            λ = $λ)
    g["expansion_coefficients"] = @benchmarkable expansion_coefficients($T_REF, $λ)
    g["scattering_matrix"] = @benchmarkable scattering_matrix($T_REF, $λ, $θs)
    g["extinction_cross_section"] = @benchmarkable extinction_cross_section($T_REF, $λ)
    g["scattering_cross_section"] = @benchmarkable scattering_cross_section($T_REF, $λ)
    g["orientation_average_analytic"] = @benchmarkable RandomOrientationTransitionMatrix($T_REF)
    g["orientation_average_numerical"] = @benchmarkable orientation_average($T_REF,
                                                                           (α, β, γ) -> 1 / (8π^2);
                                                                           Nα = 20, Nβ = 20,
                                                                           Nγ = 1) seconds=30
end

# --------------------------------------------------------------------------
# Linearization — analytical EBCM Jacobian. This is the direct baseline for
# the planned `inv(𝐐)` → `lu`/factorization-reuse optimization.
# --------------------------------------------------------------------------
let g = SUITE["linearization"] = BenchmarkGroup(["jacobian"])
    x₀ = [2.0, 1.0, 1.5, 0.02, 2π]
    vars = (:a, :c, :mᵣ, :mᵢ, :λ)
    config = (; nₘₐₓ = 6, Ng = 24)
    problem = LinearizationProblem(x₀; variables = vars) do x
        (; shape = Spheroid(x[1], x[2], complex(x[3], x[4])), λ = x[5])
    end
    g["ebcm_analytic_jacobian"] = @benchmarkable linearize_transition_matrix($problem,
                                                                             EBCMLinearization();
                                                                             config = $config) seconds=30
end

# --------------------------------------------------------------------------
# Precision — same EBCM block in Float64 vs Double64, to quantify the cost of
# the precision-escalation path the roadmap proposes to automate.
# --------------------------------------------------------------------------
let g = SUITE["precision"] = BenchmarkGroup(["multiprecision"])
    spheroid_d64 = Spheroid(Double64(2.0), Double64(1.0),
                            Complex{Double64}(Double64(1.5), Double64(0.02)))
    g["ebcm_m0_Float64"] = @benchmarkable transition_matrix_m₀($SPHEROID, $λ, 8, 32)
    g["ebcm_m0_Double64"] = @benchmarkable transition_matrix_m₀($spheroid_d64,
                                                                2 * Double64(π), 8, 32) seconds=30
end

# --------------------------------------------------------------------------
# Reuse tuned parameters when available so CI/repeat runs stay comparable.
# --------------------------------------------------------------------------
let paramspath = joinpath(@__DIR__, "params.json")
    if isfile(paramspath)
        loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals, :samples)
    end
end

SUITE
