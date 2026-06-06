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

### Support at a glance

| Backend | Model / shapes | Canonical variables | Notes |
| --- | --- | --- | --- |
| `MieLinearization()` | sphere (direct size-parameter model) | `:x`, `:mᵣ`, `:mᵢ`, `:λ` | `linearize_observable` covers `scattering_cross_section`, `extinction_cross_section`, `absorption_cross_section`, `albedo`, and fixed-angle `amplitude_matrix` |
| `EBCMLinearization()` | spheroid, Chebyshev, cylinder (fixed `nₘₐₓ`, `Ng`) | spheroid `:a`, `:c`, `:mᵣ`, `:mᵢ`, `:λ`; Chebyshev `:r₀`, `:ε`, `:mᵣ`, `:mᵢ`, `:λ`; cylinder `:r`, `:h`, `:mᵣ`, `:mᵢ`, `:λ` | forward tangents propagated through `T = -P * inv(P + im * U)` |
| `IITMLinearization()` (`:axisymmetric` / `:nfold` / `:arbitrary`) | fixed geometry | `:mᵣ`, `:mᵢ`, `:λ` | geometry variables not yet supported analytically |

Always confirm a specific combination at runtime with
`supports_linearization(problem, backend; output = :transition_matrix)`.

The Mie backend is the first implemented backend. It supports the direct size
parameter model with unique canonical variables drawn from `:x`, `:mᵣ`,
`:mᵢ`, and `:λ`.
`linearize_transition_matrix(problem, MieLinearization())` returns a
`MieTransitionMatrix` value and one derivative `MieTransitionMatrix` per
parameter. `linearize_observable` currently supports
`scattering_cross_section`, `extinction_cross_section`,
`absorption_cross_section`, `albedo`, and fixed-angle `amplitude_matrix` for
this backend. The amplitude matrix observable requires angle configuration,
for example `config = (; angles = (ϑᵢ, φᵢ, ϑₛ, φₛ))`.

The EBCM backend currently has fixed-order analytical slices for spheroids,
Chebyshev particles, and cylinders. Spheroids support unique
canonical variables drawn from `:a`, `:c`, `:mᵣ`, `:mᵢ`, and `:λ`;
Chebyshev particles support `:r₀`, `:ε`, `:mᵣ`, `:mᵢ`, and `:λ`, with the
Chebyshev degree `n` treated as a fixed discrete parameter. Cylinders support
`:r`, `:h`, `:mᵣ`, `:mᵢ`, and `:λ`; their geometry slice also differentiates
the moving piecewise quadrature nodes, weights, and Wigner angular terms.
These slices use canonical rebuild/config fields `shape`, `λ`, `nₘₐₓ`, and
`Ng`. They compute analytical forward tangents for the quadrature radius terms,
material index,
wavelength-dependent wavenumber, and Ricatti-Bessel functions, then propagate
`dP` and `dU` through the analytical identity `T = -P * inv(P + im * U)`.
`ForwardDiff.jl` is still used in tests as the reference, but not in the EBCM
production linearization path.

A spheroid EBCM Jacobian is requested the same way as the Mie one — the rebuild
function returns the canonical `shape`, `λ`, `nₘₐₓ`, and `Ng` fields:

```julia
problem = LinearizationProblem(
    [2.0, 3.0, 1.311, 0.02, 2π];
    variables = (:a, :c, :mᵣ, :mᵢ, :λ),
) do x
    (; shape = Spheroid(x[1], x[2], complex(x[3], x[4])), λ = x[5], nₘₐₓ = 10, Ng = 100)
end

result = linearize_transition_matrix(problem, EBCMLinearization())
∂T_∂a = derivative(result, :a)   # ∂T / ∂a  (semi-major axis)
```

The IITM backend currently has fixed-geometry analytical slices for
axisymmetric, n-fold, and arbitrary-shape solvers. `IITMLinearization()`,
`IITMLinearization(:axisymmetric)`, `IITMLinearization(:nfold)`, and
`IITMLinearization(:arbitrary)` support unique canonical variables drawn from
`:mᵣ`, `:mᵢ`, and `:λ`. The axisymmetric slice requires canonical
rebuild/config fields `shape`, `λ`, `nₘₐₓ`, `Nr`, and `Nϑ`; n-fold and
arbitrary-shape slices also require `Nφ`. These slices differentiate the Mie
initialization, radial Ricatti-Bessel blocks, shell interaction matrix,
projected `Q` blocks, derivative Fourier coefficients where the value solver
uses Fourier acceleration, and the IITM T-matrix recurrence. Geometry variables
remain unsupported analytically until boundary-aware shape derivatives are
implemented.

## Intended implementation layers

The framework separates the work into two layers:

- Solver kernels compute `T` and `dT/dx`. Mie establishes the first vertical
  slice. EBCM now has analytical spheroid and Chebyshev slices for fixed solver
  order and quadrature, plus a cylinder slice that includes moving quadrature
  and angular-function derivatives. IITM now has analytical fixed-geometry
  slices for axisymmetric, n-fold, and arbitrary-shape material/wavelength
  variables; the same backend API still reserves boundary-aware shape
  derivatives for later slices.
- Observable kernels propagate `T` and `dT/dx` to cross sections, albedo,
  asymmetry parameter, amplitude matrices, and scattering matrices.

For EBCM, the central matrix identity is already visible in the implementation:
`T = -P * inv(P + im * U)`. Future analytical kernels should extend the current
smooth-shape derivative pattern beyond spheroids and Chebyshev particles.
Prisms and IITM shape derivatives need separate treatment because their
indicator functions or piecewise boundaries introduce non-smooth parameter
dependence.

## Optional follow-ups

- TODO: Consider an adjoint/reverse-recurrence path for scalar Mie observables
  such as `scattering_cross_section`, `extinction_cross_section`,
  `absorption_cross_section`, `albedo`, and `asymmetry_parameter`. This path
  should avoid materializing the full transition-matrix Jacobian when users only
  need gradients of scalar observables, and it should be validated against both
  the current forward analytical linearization and `ForwardDiff.jl`.
- TODO: Extend IITM analytical linearization beyond the current fixed-geometry
  material/wavelength slices once the boundary-aware shape derivative design is
  validated.

## References

- Spurr, Wang, Zeng, and Mishchenko, "Linearized T-matrix and Mie scattering
  computations", 2012.
- Feng Xu and Anthony B. Davis, "Derivatives of light scattering properties of a
  nonspherical particle computed with the T-matrix method", Optics Letters 36(22),
  4464--4466, 2011, doi:10.1364/OL.36.004464.
- Sun, Gao, Bi, and Spurr, "Analytical Jacobians of single scattering optical
  properties using the invariant imbedding T-matrix method", Optics Express
  29(6), 9635--9669, 2021, doi:10.1364/OE.421886.
- Gao and Sun, "Improvement and application of linearized invariant imbedding
  T-matrix scattering method", Journal of Quantitative Spectroscopy and
  Radiative Transfer 290, 108322, 2022, doi:10.1016/j.jqsrt.2022.108322.
