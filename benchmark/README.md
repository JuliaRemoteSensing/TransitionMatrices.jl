# Benchmarks

Performance benchmark suite for `TransitionMatrices.jl`, built on
[BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl). The suite
follows the standard Julia convention — `benchmark/benchmarks.jl` defines a
single global `const SUITE::BenchmarkGroup` — so it works unchanged with both
[PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl) and
[AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl).

## Layout

| File | Purpose |
| --- | --- |
| `benchmarks.jl` | Defines `SUITE` (consumed by every tool below). |
| `runbenchmarks.jl` | Local runner: tune once, cache params, save results JSON. |
| `Project.toml` | Isolated environment; dev-depends on the parent package. |
| `params.json` | Cached `tune!` parameters (created on first run). |

## Groups

- `special_functions` — Riccati-Bessel and Wigner-d recursions (inner kernels).
- `ebcm` — EBCM blocks, fixed-order solves, and end-to-end auto-convergence.
- `iitm` — Invariant Imbedding T-Matrix (axisymmetric and N-fold).
- `postprocessing` — far-field observables from a precomputed T-matrix.
- `linearization` — analytical EBCM Jacobian (baseline for the `inv` → `lu` work).
- `precision` — Float64 vs Double64 on the same EBCM block.

## Running

The simplest path uses the parent `justfile`:

```sh
just bench            # instantiate the env and run the whole suite
```

Or directly:

```sh
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark benchmark/runbenchmarks.jl results-baseline.json
```

Run a single group while iterating:

```julia
julia --project=benchmark
julia> include("benchmark/benchmarks.jl");
julia> run(SUITE["ebcm"]; verbose = true)
```

## Comparing two revisions

Capture a baseline, make your change, capture again, then judge:

```julia
using BenchmarkTools
old = BenchmarkTools.load("results-baseline.json")[1]
new = BenchmarkTools.load("results.json")[1]
judge(minimum(new), minimum(old))   # :improvement / :regression / :invariant
```

`minimum` is used because it is the most noise-robust statistic for these
deterministic numerical kernels.

For automated cross-revision runs (and PR comments), AirspeedVelocity.jl reads
the same `SUITE`:

```sh
benchpkg TransitionMatrices --rev=main,HEAD --bench-on=HEAD
```
