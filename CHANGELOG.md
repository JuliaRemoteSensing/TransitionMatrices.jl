# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

No unreleased changes yet.

## [0.5.0] - 2026-06-05

### Added

- **Sh-matrix moment-separation method** for fast parameter sweeps:
  `prepare_sh`, `transition_matrix(prep, λ, mᵣ)`, and `transition_matrix_spectrum`.
  The geometry quadrature is computed once and reused, so wavelength /
  refractive-index sweeps are cheap; for spheroids it also reproduces the
  high-aspect stabilization analytically. Works for any axisymmetric shape.
- **Two-layer solver API**: fixed-discretization solvers `EBCM`, `IITM`,
  `ShMatrix`, and the `Iterative` convergence wrapper (with `ConvergencePolicy`),
  selected via `transition_matrix(s, λ, solver)`.
- **Stabilized high-aspect-ratio spheroids** (`stable = true`): the
  cancellation-free `F⁺` formulation of Somerville, Auguié & Le Ru (2013),
  available as `EBCM(nₘₐₓ, Ng; stable = true)` and, auto-converged, as
  `Iterative(EBCM; stable = true)`.
- **Performance benchmark suite** and a PR benchmark CI workflow
  (AirspeedVelocity).
- **`examples/` folder** with six runnable Pluto notebooks — shapes gallery,
  solver landscape, angular scattering & polarization, orientation averaging,
  and spectral sensitivity, plus a rain-radar observables example adapted from
  [`@xiongyuup`](https://github.com/xiongyuup)'s
  [PR #3](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/pull/3)
  with water/ice refractive-index fits and a multi-frequency
  brightness-temperature lookup — rendered into the documentation site.
- **Optional radiative-transfer workflow example** with synthetic WRF-like
  NetCDF generation and a WRF-style dual-pol radar pipeline under
  `examples/radiative_transfer/`.

### Changed

- **(Breaking)** Solver selection moves to solver objects passed to
  `transition_matrix(s, λ, solver)`. The keyword-based auto-converge
  `transition_matrix(s, λ; threshold, …)` is replaced by
  `Iterative(EBCM; threshold, …)`. The bare `transition_matrix(s, λ)` /
  `calc_T(s, λ)` is unchanged (auto-converged classic EBCM).
- Performance: factor `𝐐` once and reuse it across Jacobian slices instead of an
  explicit `inv`; the IITM radial recursion solves `M \ X` instead of forming
  `inv(M)`; the Sh-matrix precompute is optimized (coefficient precompute, folded
  convolution, precision split, and an analytic-zeroing fast path for spheroids).
- Documentation: the API reference now filters internal (underscore-prefixed)
  docstrings and uses `checkdocs = :exported`.

### Fixed

- `prepare_sh` requires an even `Ng` only for shapes with an equatorial symmetry
  plane (non-symmetric shapes accept odd `Ng`).
- `transition_matrix_spectrum` validates that `λs` and `mᵣs` have matching lengths
  instead of silently truncating to the shorter one.
- `ricattibessely` handles `nₘₐₓ = 1`.
- `cbrt` on `Double64` (DoubleFloats v1.9 returns a tuple).

## [0.4.0] - 2026-06-04

- Performance optimizations and the analytical linearization framework
  (`LinearizationProblem`, `linearize_transition_matrix`, with Mie / EBCM / IITM
  backends). See [#7](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/pull/7).

## [0.3.1] - 2026-01-18

- Raise the `Wigxjpf` compatibility bound to v0.2; CI configuration updates.

## [0.3.0] - 2025-04-03

- Compatible with Julia 1.13.
- Removed `ArbNumerics` from the dependencies; its types are no longer
  re-exported (use `Arblib`'s `Arb`/`Acb`).

## [0.2.0] - 2023-04-02

- Early development release.

## [0.1.0] - 2023-03-17

- Initial release.

[Unreleased]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/releases/tag/v0.1.0
