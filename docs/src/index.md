```@meta
CurrentModule = TransitionMatrices
```

# TransitionMatrices.jl

The transition matrix method, or T-Matrix method, is one of the most powerful and widely used tools for rigorously computing electromagnetic scattering by single and compounded particles. As a package focusing on this method, `TransitionMatrices.jl` provides the following features:

- Calculate the T-Matrix of various types of scatterers
  - Homogeneous spheres (via `bhmie`)
  - Coated spheres (via `bhcoat`)
  - Homogeneous axisymmetric shapes:
    - Spheroids
    - Cylinders
    - Chebyshev particles
- Calculate far-field scattering properties using the T-Matrix
  - Cross sections and single scattering albedo (SSA)
  - Amplitude scattering matrix
  - Phase matrix
  - Scattering matrix
