# Performance & benchmarks

## Performance notes

A few rules of thumb when scattering calculations get expensive:

- **Cost grows with the size parameter.** The truncation order needed for
  convergence scales roughly with ``x = 2\pi r/\lambda``, and the dense linear
  algebra is super-linear in that order. Large particles are intrinsically
  costly.
- **Reach for `stable` before raising precision on high-aspect spheroids.** For
  an elongated spheroid in `Float64`, `Iterative(EBCM; stable = true)` is both
  *more accurate and faster* than bumping the element type to `Double64` without
  it; see [Usage](@ref) → *High-aspect-ratio spheroids*.
- **Extended precision is several times slower.** `Double64`/`Arb` buy headroom
  against round-off and larger size parameters, but at a real runtime cost. Use
  them when accuracy demands it, not by default.
- **Reuse a preparation for parameter sweeps.** When you need one shape at many
  wavelengths or refractive indices, [`prepare_sh`](@ref) / [`ShMatrix`](@ref)
  amortise the geometry quadrature across all points.
- **Multithreading.** Several kernels use threads; start Julia with
  `julia --threads=auto` (or set `JULIA_NUM_THREADS`) to use them.

## The benchmark suite

The repository ships a benchmark suite built on
[BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl), in
[`benchmark/`](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/tree/main/benchmark).
It follows the standard Julia convention — `benchmark/benchmarks.jl` defines a
single global `const SUITE::BenchmarkGroup` — so it works unchanged with both
[PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl) and
[AirspeedVelocity.jl](https://github.com/MilesCranmer/AirspeedVelocity.jl).

Benchmark groups:

- `special_functions` — Riccati–Bessel and Wigner-d recursions (inner kernels).
- `ebcm` — EBCM blocks, fixed-order solves, and end-to-end auto-convergence.
- `iitm` — invariant-imbedding T-matrix (axisymmetric and N-fold).
- `postprocessing` — far-field observables from a precomputed T-matrix.
- `linearization` — the analytical EBCM Jacobian.
- `precision` — `Float64` vs `Double64` on the same EBCM block.

### Running

The simplest path uses the repository's `justfile`:

```sh
just bench            # instantiate the env and run the whole suite
```

Or directly, in the isolated `benchmark/` environment:

```sh
julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmark benchmark/runbenchmarks.jl results-baseline.json
```

Absolute timings are machine-specific, so no numbers are reproduced here; the
suite is meant for tracking *relative* changes (regressions and speed-ups)
between revisions. See [`benchmark/README.md`](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/blob/main/benchmark/README.md)
for the full layout.
