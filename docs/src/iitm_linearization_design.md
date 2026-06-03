# IITM Analytical Linearization Design

This note sketches the intended implementation of analytical Jacobians for the
invariant imbedding T-matrix method (IITM). Fixed-geometry material/wavelength
slices are implemented for the axisymmetric, n-fold, and arbitrary-shape
solvers; geometry derivatives remain future work.

The design follows the current package policy: analytical backends must either
compute analytical derivatives or report `UnsupportedLinearization`. They must
not hide a `ForwardDiff.jl` or finite-difference fallback in production code.

## Literature Basis

The current IITM implementation follows a shell-by-shell recurrence. The
nonlinear value path is closest to the invariant-imbedding formulation of
Johnson and to the efficient IITM/SOV implementation used by Bi, Yang,
Kattawar, and Mishchenko. Recent linearized IITM work gives the derivative
target for this package: analytical Jacobians of the T-matrix and derived
single-scattering properties with respect to refractive index, particle size,
and aspect ratio.

The most relevant papers are:

- Johnson, "Invariant imbedding T matrix approach to electromagnetic
  scattering", Applied Optics 27(23), 4861--4873, 1988,
  doi:10.1364/AO.27.004861.
- Bi, Yang, Kattawar, and Mishchenko, "Efficient implementation of the
  invariant imbedding T-matrix method and the separation of variables method
  applied to large nonspherical inhomogeneous particles", Journal of
  Quantitative Spectroscopy and Radiative Transfer 116, 169--183, 2013,
  doi:10.1016/j.jqsrt.2012.11.014.
- Sun, Gao, Bi, and Spurr, "Analytical Jacobians of single scattering optical
  properties using the invariant imbedding T-matrix method", Optics Express
  29(6), 9635--9669, 2021, doi:10.1364/OE.421886.
- Gao and Sun, "Improvement and application of linearized invariant imbedding
  T-matrix scattering method", Journal of Quantitative Spectroscopy and
  Radiative Transfer 290, 108322, 2022, doi:10.1016/j.jqsrt.2022.108322.
- Sun, Gao, Liang, Liu, and Liu, "Capability and convergence of linearized
  invariant-imbedding T-matrix and physical-geometric optics methods for light
  scattering", Optics Express 30(21), 37769--37785, 2022,
  doi:10.1364/OE.473075.

The 2021 paper is the primary target for linearized IITM derivatives. The 2022
JQSRT paper matters because it discusses improved derivatives for effective
radius definitions and aspect-ratio sensitivity, including the divergent cases
that a package implementation should avoid by construction.

## Current Code Boundary

The package currently has three IITM value solvers:

- `transition_matrix_iitm(::AbstractAxisymmetricShape, λ, nₘₐₓ, Nr, Nϑ)`
- `transition_matrix_iitm(::AbstractNFoldShape, λ, nₘₐₓ, Nr, Nϑ, Nφ)`
- `transition_matrix_iitm(::AbstractShape, λ, nₘₐₓ, Nr, Nϑ, Nφ)`

All three use the same logical recurrence:

1. Build radial quadrature between `rmin(shape)` and `rmax(shape)`.
2. Initialize the embedded T-matrix from Mie coefficients of the inscribed
   sphere at `rmin(shape)`.
3. For each radial shell, build `𝐉`, `𝐇`, `𝐆`, and `𝐔`.
4. Build the shell interaction matrix `𝐐`.
5. Update the accumulated T-matrix.
6. Repack the internal block layout into `AxisymmetricTransitionMatrix` or
   `TransitionMatrix`.

The analytical implementation should therefore live in a new
`src/IITM/linearization.jl` included from `src/IITM/index.jl`. It should not
fork a separate public model of IITM. Instead, the value solver and the
linearized solver should share private kernels for quadrature, Bessel tables,
Wigner tables, Fourier coefficients, shell matrices, and output packing.

## Public Contract

The public entry point remains:

```julia
linearize_transition_matrix(problem, IITMLinearization(variant); config)
```

where `variant` is one of `:axisymmetric`, `:nfold`, `:arbitrary`, or `:auto`.
The `config` or rebuilt problem must provide:

- `shape`
- `λ`
- `nₘₐₓ`
- `Nr`
- `Nϑ`
- `Nφ` for `:nfold` and `:arbitrary`

The first implementation should support unique canonical variables drawn from
the shape-specific set. It should not add aliases such as `:size`,
`:aspect_ratio`, `:effective_radius`, or `:shape_factor`.

Planned canonical sets:

- `Spheroid`: `:a`, `:c`, `:mᵣ`, `:mᵢ`, `:λ`
- `Cylinder`: `:r`, `:h`, `:mᵣ`, `:mᵢ`, `:λ`
- `Chebyshev`: no first IITM geometry slice until a volume-domain derivative
  is designed for its radius perturbation; material and wavelength derivatives
  can still be considered if the value path is valid.
- `Prism`: `:a`, `:h`, `:mᵣ`, `:mᵢ`, `:λ`
- General inhomogeneous `AbstractShape`: initially `:mᵣ`, `:mᵢ`, and `:λ`
  only, unless the shape provides an explicit derivative provider.

The variable tuple may contain any subset of the supported canonical variables,
in any order. Solver order and quadrature counts are fixed discrete settings.
They are not differentiated.

Derived parameters from the literature, such as aspect ratio or the effective
radii `r_V`, `r_S`, and `r_VS`, should be handled by an explicit chain rule over
canonical derivatives, not by adding multiple backend aliases for the same
physical parameter. A later helper can expose those chain-rule combinations if
that becomes common enough.

## Linearized Recurrence

For each radial shell, the value solver computes:

```text
𝐐 = wᵣ * inv(𝐈 - wᵣ * 𝐔 * 𝐆) * 𝐔
𝐐ⱼⱼ = im * k * transpose(𝐉) * 𝐐 * 𝐉
𝐐ⱼₕ = im * k * transpose(𝐉) * 𝐐 * 𝐇
𝐐ₕⱼ = im * k * transpose(𝐇) * 𝐐 * 𝐉
𝐐ₕₕ = im * k * transpose(𝐇) * 𝐐 * 𝐇
𝐓⁺ = 𝐐ⱼⱼ + (𝐈 + 𝐐ⱼₕ) * inv(𝐈 - 𝐓 * 𝐐ₕₕ) * 𝐓 * (𝐈 + 𝐐ₕⱼ)
```

The tangent implementation should propagate one directional derivative `∂`
through the same recurrence. For the shell interaction matrix, define:

```text
𝐑 = 𝐈 - wᵣ * 𝐔 * 𝐆
𝐑 * ∂𝐐 = ∂wᵣ * 𝐔 + wᵣ * ∂𝐔 - ∂𝐑 * 𝐐
∂𝐑 = -∂wᵣ * 𝐔 * 𝐆 - wᵣ * ∂𝐔 * 𝐆 - wᵣ * 𝐔 * ∂𝐆
```

This form is preferred over differentiating `inv(𝐑)` directly because it can
reuse the same factorization used for the value solve.

Each projected `𝐐` block follows the product rule. For example:

```text
∂𝐐ⱼₕ =
    im * ∂k * transpose(𝐉) * 𝐐 * 𝐇 +
    im * k  * transpose(∂𝐉) * 𝐐 * 𝐇 +
    im * k  * transpose(𝐉) * ∂𝐐 * 𝐇 +
    im * k  * transpose(𝐉) * 𝐐 * ∂𝐇
```

For the T recurrence, define:

```text
𝐀 = 𝐈 + 𝐐ⱼₕ
𝐁 = 𝐈 - 𝐓 * 𝐐ₕₕ
𝐂 = 𝐈 + 𝐐ₕⱼ
𝐒 = inv(𝐁)
```

Then:

```text
∂𝐁 = -∂𝐓 * 𝐐ₕₕ - 𝐓 * ∂𝐐ₕₕ
∂𝐒 = -𝐒 * ∂𝐁 * 𝐒
∂𝐓⁺ =
    ∂𝐐ⱼⱼ +
    ∂𝐀 * 𝐒 * 𝐓 * 𝐂 +
    𝐀 * ∂𝐒 * 𝐓 * 𝐂 +
    𝐀 * 𝐒 * ∂𝐓 * 𝐂 +
    𝐀 * 𝐒 * 𝐓 * ∂𝐂
```

This is a forward analytical recurrence. It scales linearly in the number of
requested variables after the value matrices and factorizations are available.
An adjoint or reverse recurrence for scalar observables can be considered later,
but it should not block the first T-matrix Jacobian implementation.

## Local Derivative Sources

The recurrence above depends on local tangents from four sources.

### Material and Wavelength

For all variants:

```text
k = 2π / λ
∂k = -k * ∂λ / λ
∂m = 1       for :mᵣ
∂m = im      for :mᵢ
∂m = 0       otherwise
```

Inside homogeneous particles:

```text
ε = m^2
∂ε = 2m * ∂m
∂((ε - 1) / ε) = ∂ε / ε^2
```

For inhomogeneous shapes, the first implementation should require either a
constant material index or a shape-provided `refractive_index` derivative. It
should not infer derivatives through arbitrary user code.

### Radial Quadrature

Radial nodes and weights depend on `rₘᵢₙ` and `rₘₐₓ`:

```text
rᵢ = (rₘₐₓ - rₘᵢₙ) * (xᵢ + 1) / 2 + rₘᵢₙ
wᵢ = (rₘₐₓ - rₘᵢₙ) * ŵᵢ / 2
```

Their tangents are direct product-rule derivatives. Shape variables are
supported only away from branch points of `min`, `max`, and geometric partition
changes. For example, a spheroid with `a == c` has a nonsmooth `rmin/rmax`
branch and should report unsupported shape derivatives unless a dedicated
spherical path is used.

### Bessel and Hankel Blocks

`𝐉`, `𝐇`, and `𝐆` depend on `k r`. The existing EBCM analytical code already
uses the Ricatti-Bessel differential equation pattern:

```text
∂f = f′ * ∂z
∂f′ = (n(n+1)/z^2 - 1) * f * ∂z
```

The IITM linearization should reuse or move that helper into
`src/special_functions`, so EBCM and IITM do not maintain separate copies.

### Wigner and Angular Tables

The value solvers precalculate Wigner `d`, `𝜋`, and `τ` tables. If angular
nodes are fixed, their tangents are zero. If a geometry implementation uses
moving angular partitions, the same derivative identity already used by EBCM
cylinder linearization applies:

```text
∂d = τ * ∂ϑ
∂τ = d² * ∂ϑ
∂𝜋 = m * (∂d / sinϑ - d * cosϑ * ∂ϑ / sinϑ^2)
```

where:

```text
d² = -cosϑ / sinϑ * τ - (n(n+1) - m^2 / sinϑ^2) * d
```

This helper should also be shared with EBCM rather than duplicated.

## Geometry Derivatives

Geometry is the main risk. The current IITM value code evaluates pointwise
membership through `x ∈ shape` or `refractive_index(shape, x)`. That is fine for
values, but it is not a valid analytical shape derivative: the derivative of a
moving indicator is a boundary distribution, not a pointwise derivative of the
Boolean result.

The implementation should therefore split geometry support into three tiers.

### Tier 1: Fixed Geometry

Support `:mᵣ`, `:mᵢ`, and `:λ` with fixed geometry for axisymmetric, n-fold,
and arbitrary variants. This requires no derivative of membership. It exercises
the full shell recurrence and is the right first vertical slice.

### Tier 2: Uniform Size Scaling

Support a uniform size parameter only when the shape can be represented as a
scaled reference body. If `rₘᵢₙ`, `rₘₐₓ`, and radial nodes scale together, the
membership of a normalized quadrature point is invariant. The derivative then
comes from moving radial nodes, radial weights, `kr`, and the initial inscribed
sphere, not from differentiating a Boolean indicator.

This tier is useful for effective-radius sensitivities, but the backend should
still expose canonical derivatives first. Derived effective-radius gradients can
be formed by chain rule.

### Tier 3: Aspect and Boundary Shape Derivatives

Aspect-ratio derivatives require boundary-aware integration. Fixed angular
masks are not acceptable.

For axisymmetric shapes, the practical path is a
`LinearizedIITMGeometry` abstraction that can return, for each radial shell:

- volume angular quadrature nodes and weights,
- tangents of any moving angular nodes and weights,
- Wigner table tangents induced by moving angular nodes,
- material contrast coefficients and their tangents,
- optional surface-correction terms if the derivative is represented as a
  boundary integral.

Spheroids and cylinders can then be implemented with shape-specific angular
partitions. Cylinders are piecewise smooth and require separate side/cap
branches. Spheroids are smooth but still need root/partition handling when a
radial shell intersects the boundary.

For n-fold prisms, the same idea extends to moving azimuthal partitions or a
surface-quadrature correction over facets. This should be a later phase because
it touches the Fourier coefficient path and the block-decoupled n-fold solver.

For arbitrary shapes, analytical shape derivatives should require an explicit
shape derivative provider. The package should not try to infer boundary motion
from an arbitrary `∈` method.

## Fourier Coefficients

The n-fold and arbitrary solvers use Fourier coefficients of:

```text
ε - 1
(ε - 1) / ε
```

The linearized path needs matching coefficients for their tangents. In the
`ComplexF64` FFT path, this means extending `_azimuthal_fourier_coefficients`
to accept optional derivative work arrays and return:

```text
coeff_ε, coeff_εinv, ∂coeff_ε, ∂coeff_εinv
```

The non-FFT path can accumulate value and tangent in the same explicit
quadrature loops. Both paths must produce identical mathematical quantities, so
tests should cover both by using a non-`ComplexF64` element type in at least one
small case.

## Initialization Derivative

The initial embedded T-matrix comes from Mie coefficients at:

```text
x₀ = k * rₘᵢₙ
```

Its tangent is:

```text
∂x₀ = ∂k * rₘᵢₙ + k * ∂rₘᵢₙ
```

The implementation should reuse the existing Mie analytical derivative kernel
instead of re-deriving coefficient derivatives in IITM. This keeps the sphere
limit consistent across backends.

## Implementation Phases

1. Refactor IITM value kernels so value and tangent code can share quadrature,
   Bessel, Wigner, Fourier, shell matrix, and output-packing utilities. This
   phase should preserve existing numerical results.
2. Implement `IITMLinearization(:axisymmetric)` for fixed geometry variables
   `:mᵣ`, `:mᵢ`, and `:λ`.
3. Add fixed-geometry support for `:nfold` and `:arbitrary`, including
   derivative Fourier coefficients.
4. Add uniform size scaling for shapes with a well-defined scaling derivative.
   Validate the sphere limit against Mie.
5. Add axisymmetric boundary-aware geometry derivatives for `Spheroid` and
   `Cylinder`, using canonical variables `:a`, `:c`, `:r`, and `:h`.
6. Add n-fold prism geometry derivatives after the Fourier coefficient and
   moving-boundary design is validated on axisymmetric shapes.
7. Consider adjoint/reverse recurrences for scalar observables after the
   forward T-matrix Jacobian is correct.

## Verification Plan

Tests should separate value agreement, tangent agreement, and unsupported
boundaries.

- Material and wavelength tangents: compare with `ForwardDiff.jl` where the
  value path supports dual numbers, and with central finite differences
  otherwise.
- Sphere limit: compare IITM `:r`, `:mᵣ`, `:mᵢ`, and `:λ` derivatives against
  Mie analytical derivatives.
- Axisymmetric shapes: compare spheroid and cylinder derivatives against EBCM
  analytical derivatives at small `nₘₐₓ`, moderate quadrature, and nonsingular
  geometry.
- Literature-level behavior: reproduce at least one qualitative sensitivity
  case from Sun et al. 2021 or Gao and Sun 2022 after shape derivatives land.
- Unsupported cases: assert that nonsmooth branch points, unknown shape
  derivative providers, duplicated variables, and unsupported aliases fail
  explicitly.

For shape derivatives, `ForwardDiff.jl` is not a sufficient reference if the
value path contains Boolean membership tests. Finite differences and cross
backend comparisons are required.

## Performance Expectations

The value IITM cost is dominated by shell-local matrix construction,
factorization, and T recurrence updates. A forward analytical Jacobian adds one
tangent recurrence per requested variable, but it can share:

- radial and angular quadrature,
- Wigner value tables,
- Bessel value tables,
- FFT plans and, where possible, value spectra,
- shell matrix factorizations for `𝐑` and `𝐁`.

The expected complexity is approximately:

```text
cost ≈ value_cost + nvars * tangent_cost
```

This is still forward-mode complexity, but it avoids the generic dual-number
overhead and accidental allocation growth of `ForwardDiff` through large
complex matrices. For scalar observables with many variables, a later adjoint
recurrence could avoid materializing the full T-matrix Jacobian.

## Open Risks

- Shape derivatives through fixed pointwise membership are mathematically
  wrong. Geometry support must be boundary-aware.
- `rmin` and `rmax` are nonsmooth at branch changes such as `a == c` or
  `r == h/2`.
- Existing IITM code uses `inv` in several places. The linearized path should
  prefer factorization/solve forms, but changing the value path must be tested
  separately.
- Threaded loops need per-thread tangent workspaces. Shared derivative buffers
  would introduce races.
- The `ComplexF64` FFT fast path and generic explicit quadrature path must keep
  the same coefficient normalization.
- Shape variables from the literature are often effective radius and aspect
  ratio, while this package should keep canonical geometric variables in the
  backend. The chain rule layer must be explicit to avoid alias complexity.
