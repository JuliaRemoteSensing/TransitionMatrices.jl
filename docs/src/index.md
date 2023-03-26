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

Compared to existing packages, `TransitionMatrices.jl` is special in that it is generic and supports various floating-point types, e.g.:

- `Float64` and `BigFloat` from [`Base`](https://docs.julialang.org/en/v1/base/)
- `Double64` from [`DoubleFloats.jl`](https://github.com/JuliaMath/DoubleFloats.jl)
- `Float128` from [`Quadmath.jl`](https://github.com/JuliaMath/Quadmath.jl)
- `Arb` from [`Arblib.jl`](https://github.com/kalmarek/Arblib.jl)
- `ArbFloat` from [`ArbNumerics.jl`](https://github.com/JeffreySarnoff/ArbNumerics.jl)

For EBCM, by using higher-precision floating-point types, the maximum size parameter that can be handled is greatly improved.
