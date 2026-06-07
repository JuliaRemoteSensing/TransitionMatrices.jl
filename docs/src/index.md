```@meta
CurrentModule = TransitionMatrices
```

# TransitionMatrices.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaRemoteSensing.github.io/TransitionMatrices.jl/dev/)
[![Build Status](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaRemoteSensing/TransitionMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaRemoteSensing/TransitionMatrices.jl)

The transition matrix method, or **T-Matrix method**, is one of the most powerful
and widely used tools for rigorously computing electromagnetic scattering by
single and compounded particles. `TransitionMatrices.jl` is a generic,
arbitrary-precision Julia implementation focused on this method.

## Installation

`TransitionMatrices.jl` requires Julia ≥ 1.10. From the Julia REPL's package mode
(press `]`):

```julia-repl
pkg> add TransitionMatrices
```

To track the development version, or if the package is not yet in your registry,
add it by URL — this always works:

```julia-repl
pkg> add https://github.com/JuliaRemoteSensing/TransitionMatrices.jl
```

## Quick start

```julia
using TransitionMatrices

# A prolate spheroid: semi-axes a = 1, c = 2, relative refractive index m = 1.5 + 0.01im
s = Spheroid(1.0, 2.0, 1.5 + 0.01im)

# Auto-converged (classic EBCM) T-matrix at wavelength λ = 2π
𝐓 = transition_matrix(s, 2π)

# Orientation-averaged far-field observables
Qsca = scattering_cross_section(𝐓, 2π)   # scattering cross section
Qext = extinction_cross_section(𝐓, 2π)   # extinction cross section
ω    = albedo(𝐓)                         # single-scattering albedo
g    = asymmetry_parameter(𝐓, 2π)        # asymmetry parameter ⟨cos Θ⟩
```

That is the whole loop: **define a shape → build its T-matrix → derive far-field
quantities**. From here, [Usage](@ref) walks through every step, [Theory &
conventions](@ref) defines the quantities and conventions used, and [Examples](@ref)
collects runnable notebooks.

## Features

- Calculate the T-Matrix of various types of scatterers
  - Homogeneous spheres (via `bhmie`)
  - Coated spheres (via `bhcoat`)
  - Homogeneous axisymmetric shapes (via EBCM and IITM)
    - Spheroids
    - Cylinders
    - Chebyshev particles
  - Arbitrary shapes (via IITM)
    - Prisms
- Calculate far-field scattering properties using the T-Matrix
  - Cross sections and single scattering albedo (SSA)
  - Amplitude scattering matrix
  - Phase matrix
  - Scattering matrix
- Compute Jacobians through the linearization framework
  - Numerical automatic differentiation for user-defined scalar workflows via `ForwardDiff.jl`
  - Analytical Mie linearization for size, refractive-index, and wavelength variables
  - Analytical EBCM slices for spheroids, cylinders, and Chebyshev particles
  - Analytical fixed-geometry IITM material/wavelength slices for axisymmetric,
    n-fold, and arbitrary-shape solvers

## Floating-point genericity

Compared to existing packages, `TransitionMatrices.jl` is special in that it is
generic and supports various floating-point types, e.g.:

- `Float64` and `BigFloat` from [`Base`](https://docs.julialang.org/en/v1/base/)
- `Double64` from [`DoubleFloats.jl`](https://github.com/JuliaMath/DoubleFloats.jl) when `DoubleFloats` is loaded
- `Float128` from [`Quadmath.jl`](https://github.com/JuliaMath/Quadmath.jl) when `Quadmath` is loaded
- `Arb` and `Acb` from [`Arblib.jl`](https://github.com/kalmarek/Arblib.jl)

For EBCM, using a higher-precision floating-point type greatly improves the
maximum size parameter that can be handled (see [Choosing a solver](@ref) for
when to reach for extended precision versus the `stable` formulation).

The precision types `Arb` and `Acb` are re-exported by `TransitionMatrices.jl`.
The `DoubleFloats.jl` and `Quadmath.jl` precision types are optional: load them
explicitly with `using DoubleFloats` or `using Quadmath` before constructing
`Double64`, `Float128`, or `ComplexF128` values.

The `0.6` compatibility line keeps `Quadmath.jl` 1.x and `Wigxjpf.jl` 0.3.x,
and moves `DoubleFloats.jl`, `Quadmath.jl`, and `GenericFFT.jl` to optional
package extensions.

## Where to go next

- [Theory & conventions](@ref) — the T-matrix definition, the index/sign
  conventions, the size parameter, and how each observable is defined.
- [Usage](@ref) — define a shape, choose a solver, build the T-matrix, and
  compute far-field quantities, orientation averages, and derivatives.
- [Examples](@ref) — runnable [Pluto.jl](https://plutojl.org) notebooks with plots.
- [Linearization Framework](@ref) — analytical and automatic differentiation.
- [API](@ref) — the full reference, grouped by topic.
- [Methods & references](@ref) — the literature each method is built on.

## How to cite

If you use `TransitionMatrices.jl` in your research, please cite it — a
[`CITATION.bib`](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/blob/main/CITATION.bib)
is provided in the repository. Please also cite the original publication(s) for
the specific method you use; see [Methods & references](@ref).
