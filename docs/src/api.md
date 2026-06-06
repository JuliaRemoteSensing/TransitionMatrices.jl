```@meta
CurrentModule = TransitionMatrices
```

# API

The full public reference, grouped by topic. Internal (`_`-prefixed) helpers are
omitted. Use the index below to jump to any symbol.

```@index
```

## Shapes

Scatterer geometries and their common interface.

```@autodocs
Modules = [TransitionMatrices]
Pages = ["shapes/index.jl", "shapes/spheroid.jl", "shapes/cylinder.jl",
    "shapes/chebyshev.jl", "shapes/prism.jl"]
Filter = t -> !startswith(string(Base.nameof(t)), "_")
```

## Solvers & T-matrix construction

The two-layer solver API ([`AbstractSolver`](@ref), the fixed solvers, and
[`Iterative`](@ref)) and the underlying engines (Mie, EBCM with its `stable` and
Sh-matrix variants, IITM).

```@autodocs
Modules = [TransitionMatrices]
Pages = ["solvers.jl",
    "EBCM/index.jl", "EBCM/axisymmetric.jl", "EBCM/routines.jl",
    "EBCM/shmatrix.jl", "EBCM/stabilization.jl",
    "IITM/index.jl", "IITM/axisymmetric.jl", "IITM/arbitrary.jl",
    "IITM/nfold.jl", "IITM/fourier.jl",
    "Mie/index.jl", "Mie/bhmie.jl", "Mie/bhcoat.jl", "Mie/MieTransitionMatrix.jl"]
Filter = t -> !startswith(string(Base.nameof(t)), "_")
```

## Transition matrices, observables & post-processing

The transition-matrix types and the far-field quantities derived from them
(cross sections, amplitude/phase/scattering matrices, orientation averaging).

```@autodocs
Modules = [TransitionMatrices]
Pages = ["common/index.jl", "common/AbstractTransitionMatrix.jl",
    "common/AxisymmetricTransitionMatrix.jl",
    "common/RandomOrientationTransitionMatrix.jl", "common/utils.jl"]
Filter = t -> !startswith(string(Base.nameof(t)), "_")
```

## Near-field reconstruction

Vector spherical wave functions and the external electromagnetic field
(incident, scattered, total) reconstructed from any transition matrix (see the
[Near-field maps from a T-matrix](examples/near_field.md) example).

```@autodocs
Modules = [TransitionMatrices]
Pages = ["common/vswf.jl", "common/near_field.jl"]
Filter = t -> !startswith(string(Base.nameof(t)), "_")
```

## Linearization

The differentiation framework and its backends (see
[Linearization Framework](@ref) for the narrative).

```@autodocs
Modules = [TransitionMatrices]
Pages = ["linearization.jl", "Mie/linearization.jl",
    "EBCM/linearization.jl", "IITM/linearization.jl"]
Filter = t -> !startswith(string(Base.nameof(t)), "_")
```

## Special functions

The Riccati–Bessel, Wigner-d/D, Clebsch–Gordan, factorial, and quadrature
kernels used internally (a few are exported for convenience).

```@autodocs
Modules = [TransitionMatrices]
Pages = ["special_functions/index.jl", "special_functions/bessel.jl",
    "special_functions/clebschgordan.jl", "special_functions/factorial.jl",
    "special_functions/quadrature.jl", "special_functions/wignerd.jl",
    "compat/index.jl"]
Filter = t -> !startswith(string(Base.nameof(t)), "_")
```
