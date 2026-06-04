# TransitionMatrices.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaRemoteSensing.github.io/TransitionMatrices.jl/dev/)
[![Build Status](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaRemoteSensing/TransitionMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaRemoteSensing/TransitionMatrices.jl)

The transition matrix method, or T-Matrix method, is one of the most powerful and widely used tools for rigorously computing electromagnetic scattering by single and compounded particles. As a package focusing on this method, `TransitionMatrices.jl` provides the following features:

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

Compared to existing packages, `TransitionMatrices.jl` is special in that it is generic and supports various floating-point types, e.g.:

- `Float64` and `BigFloat` from [`Base`](https://docs.julialang.org/en/v1/base/)
- `Double64` from [`DoubleFloats.jl`](https://github.com/JuliaMath/DoubleFloats.jl)
- `Float128` from [`Quadmath.jl`](https://github.com/JuliaMath/Quadmath.jl)
- `Arb` and `Acb` from [`Arblib.jl`](https://github.com/kalmarek/Arblib.jl)

By using higher-precision floating-point types, the maximum size parameter that can be handled is greatly improved.

The precision types `Double64`, `Float128`, `ComplexF128`, `Arb`, and `Acb`
are re-exported by `TransitionMatrices.jl` and can be directly used after
`using TransitionMatrices`.

The `0.4` compatibility line uses `Quadmath.jl` 1.x and `Wigxjpf.jl` 0.3.x.
