# Linearization Framework

`TransitionMatrices.jl` currently supports numerical automatic differentiation
through `ForwardDiff.jl`. The linearization framework provides a separate
foundation for analytical Jacobians, following the terminology used by
linearized T-matrix literature.

The first version gives the whole package one way to describe differentiable
scattering problems, one way to return values and Jacobians, and one explicit
policy for unsupported analytical backends. It does not silently fall back to
finite differences.

## Problem model

Analytical backends operate on a `LinearizationProblem`. The problem stores a
real parameter vector and a rebuild function. This keeps the same mental model
as existing numerical differentiation examples:

```julia
problem = LinearizationProblem(
    [1.7, 1.311, 0.02, 2π];
    variables = (:x, :mᵣ, :mᵢ, :λ),
) do x
    (; x = x[1], m = complex(x[2], x[3]), λ = x[4], nₘₐₓ = 5)
end
```

The rebuild function is deliberately user-provided. Shape parameters, material
parameters, wavelength, and solver settings are not forced into a new shape
type. A backend can inspect the rebuilt input and decide whether an analytical
Jacobian is available.

## Result model

Analytical calculations return `LinearizationResult(value, jacobian, variables)`.
The `value` can be a transition matrix, a scalar cross section, an amplitude
matrix, or a scattering matrix. The `jacobian` is with respect to the real
parameter vector in the problem.

For array outputs, the preferred storage convention is to place the parameter
axis last. For transition-matrix-like outputs, implementations may also store a
vector of derivative matrices. Use `derivative(result, i)` or
`derivative(result, :λ)` to access one parameter derivative without depending on
the concrete Jacobian storage.

## Backends and support

The backend markers are:

- `MieLinearization()`
- `EBCMLinearization()`
- `IITMLinearization()`, or `IITMLinearization(:axisymmetric)`,
  `IITMLinearization(:nfold)`, `IITMLinearization(:arbitrary)`

Call `supports_linearization(problem, backend; output=:transition_matrix)` to
query whether an analytical implementation exists. If the requested combination
is not implemented, `linearize_transition_matrix` and `linearize_observable`
throw `UnsupportedLinearization`.

This explicit failure is part of the design. It prevents users from mistaking a
numerical finite-difference fallback for an analytical Jacobian.

The Mie backend is the first implemented backend. It supports the direct size
parameter model with canonical variables `:x`, `:mᵣ`, `:mᵢ`, and `:λ`.
`linearize_transition_matrix(problem, MieLinearization())` returns a
`MieTransitionMatrix` value and one derivative `MieTransitionMatrix` per
parameter. `linearize_observable` currently supports
`scattering_cross_section`, `extinction_cross_section`,
`absorption_cross_section`, `albedo`, and fixed-angle `amplitude_matrix` for
this backend. The amplitude matrix observable requires angle configuration,
for example `config = (; angles = (ϑᵢ, φᵢ, ϑₛ, φₛ))`.

## Intended implementation layers

The framework separates the work into two layers:

- Solver kernels compute `T` and `dT/dx`. Mie establishes the first vertical
  slice. EBCM is the next analytical target. IITM variants are represented by
  the same API but can remain unsupported until their linearized recursions are
  implemented.
- Observable kernels propagate `T` and `dT/dx` to cross sections, albedo,
  asymmetry parameter, amplitude matrices, and scattering matrices.

For EBCM, the central matrix identity is already visible in the implementation:
`T = -P * inv(P + im * U)`. Future analytical kernels should linearize `P` and
`U`, then propagate through this identity. Smooth axisymmetric shapes such as
spheroids and Chebyshev particles are the first practical shape targets.
Cylinders, prisms, and IITM shape derivatives need separate treatment because
their indicator functions or piecewise boundaries introduce non-smooth
parameter dependence.

## Optional follow-ups

- TODO: Consider an adjoint/reverse-recurrence path for scalar Mie observables
  such as `scattering_cross_section`, `extinction_cross_section`,
  `absorption_cross_section`, `albedo`, and `asymmetry_parameter`. This path
  should avoid materializing the full transition-matrix Jacobian when users only
  need gradients of scalar observables, and it should be validated against both
  the current forward analytical linearization and `ForwardDiff.jl`.

## References

- Spurr, Wang, Zeng, and Mishchenko, "Linearized T-matrix and Mie scattering
  computations", 2012.
- Xu and Davis, "Analytical derivatives of T-matrix light scattering
  quantities", 2011.
- Recent linearized IITM work on analytical Jacobians for particle size and
  aspect-ratio sensitivity.
