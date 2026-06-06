```@meta
CurrentModule = TransitionMatrices
```

# Theory & conventions

This page summarises the definitions and conventions used throughout
`TransitionMatrices.jl`. They follow Mishchenko, Travis & Lacis, *Scattering,
Absorption, and Emission of Light by Small Particles* (Cambridge University
Press, 2002) — referred to below as **MTL** — which the in-source docstrings
cite by equation number. Full references are listed in
[Methods & references](@ref).

## The transition matrix

For a single scatterer embedded in a non-absorbing host medium, both the
incident and the scattered electric fields are expanded in **vector spherical
wave functions** (VSWFs). The incident field uses the *regular* VSWFs
``\mathrm{Rg}\mathbf{M}_{mn}, \mathrm{Rg}\mathbf{N}_{mn}`` with expansion
coefficients ``(a_{mn}, b_{mn})``; the scattered field uses the *outgoing*
VSWFs ``\mathbf{M}_{mn}, \mathbf{N}_{mn}`` with coefficients ``(p_{mn}, q_{mn})``.
Because Maxwell's equations are linear, the two sets of coefficients are linearly
related, and the matrix of that relation is the **transition matrix** ``T`` (MTL §5):

```math
\begin{bmatrix} p_{mn} \\ q_{mn} \end{bmatrix}
  = \sum_{n'=1}^{\infty} \sum_{m'=-n'}^{n'}
    \begin{bmatrix} T^{11}_{mn m'n'} & T^{12}_{mn m'n'} \\
                    T^{21}_{mn m'n'} & T^{22}_{mn m'n'} \end{bmatrix}
    \begin{bmatrix} a_{m'n'} \\ b_{m'n'} \end{bmatrix}.
```

The transition matrix depends only on the particle (its shape, size, and
refractive index relative to the host) and on the wavelength — **not** on the
incident field. Computing it once therefore gives access to scattering for *any*
incidence direction and polarisation, and to orientation averages in closed form.

### Storage

[`AbstractTransitionMatrix`](@ref) stores ``T`` as a six-dimensional array
indexed ``T_{mn m'n'}^{kl}``: the unprimed ``(m, n)`` label the scattered
harmonic, the primed ``(m', n')`` the incident one, and ``k, l \in \{1, 2\}``
select the magnetic (``\mathbf{M}``) and electric (``\mathbf{N}``) parts. The
degree runs ``1 \le n \le n_{\max}`` and the order ``-n \le m \le n``, so a
truncation order ``n_{\max}`` fixes the matrix size. Concrete subtypes
([`TransitionMatrix`](@ref), [`AxisymmetricTransitionMatrix`](@ref),
[`MieTransitionMatrix`](@ref), [`RandomOrientationTransitionMatrix`](@ref))
exploit symmetry to store and apply ``T`` more compactly.

## Size parameter and units

Lengths and the wavelength are expressed in the same (arbitrary) unit; only their
ratio matters. The wavenumber in the host medium is ``k = 2\pi/\lambda``, and the
dimensionless **size parameter** of a scatterer of characteristic radius ``r`` is

```math
x = k r = \frac{2\pi r}{\lambda}.
```

The refractive index `m` passed to a shape is **relative** to the host medium
(``m = m_\text{particle}/m_\text{host}``); a complex `m` carries absorption in
its imaginary part. The default wavelength in the API is ``\lambda = 2\pi``, for
which ``k = 1`` and the size parameter equals the radius numerically — convenient
for quick experiments.

The truncation order needed for convergence grows with the size parameter
(roughly ``n_{\max} \gtrsim x + c\,x^{1/3}``), which is why large particles are
expensive and eventually require extended precision or the `stable` formulation;
see [Choosing a solver](@ref).

## Orientation

Because ``T`` is an operator between VSWF bases, a rotation of the particle by
Euler angles ``(\alpha, \beta, \gamma)`` acts on it through Wigner ``D``-functions
(MTL Eq. (5.29)); [`rotate`](@ref) returns the rotated transition matrix without
recomputing it from scratch.

For randomly oriented particles the orientation-averaged transition matrix has a
closed form (MTL Eq. (5.96)), implemented by
[`RandomOrientationTransitionMatrix`](@ref); [`orientation_average`](@ref) instead
integrates numerically against a user-supplied orientation distribution and
converges to the analytic result. The orientation-averaged **scattering matrix**
is expanded in generalised spherical functions with coefficients computed by
[`expansion_coefficients`](@ref).

## Far-field observables

All of the following are methods of an [`AbstractTransitionMatrix`](@ref) (see
[Post-processing](@ref) for runnable examples).

### Amplitude (Jones) matrix

For a fixed incidence direction ``\hat{\mathbf{n}}^{\text{inc}}`` and scattering
direction ``\hat{\mathbf{n}}^{\text{sca}}``, the ``2\times 2`` complex amplitude
matrix ``\mathbf{S}`` maps the incident to the scattered transverse field
(MTL Eqs. (5.11)–(5.17)):

```math
\begin{bmatrix} E_\vartheta^{\text{sca}} \\ E_\varphi^{\text{sca}} \end{bmatrix}
  = \frac{e^{ikr}}{r}\,\mathbf{S}\,
    \begin{bmatrix} E_\vartheta^{\text{inc}} \\ E_\varphi^{\text{inc}} \end{bmatrix}.
```

[`amplitude_matrix`](@ref) returns ``\mathbf{S}`` for given incidence/scattering
angles.

### Phase (Mueller) matrix

The real ``4\times 4`` phase matrix ``\mathbf{Z}`` propagates the Stokes vector
and follows from ``\mathbf{S}`` (MTL Eqs. (2.106)–(2.121)); [`phase_matrix`](@ref)
builds it from an amplitude matrix.

### Cross sections, albedo, asymmetry

The orientation-averaged scattering and extinction cross sections (MTL
Eqs. (5.140)/(5.141) and (5.102)/(5.107)) are returned by
[`scattering_cross_section`](@ref) and [`extinction_cross_section`](@ref); their
difference is the [`absorption_cross_section`](@ref). The single-scattering
[`albedo`](@ref) is ``\omega = C_\text{sca}/C_\text{ext}`` and the
[`asymmetry_parameter`](@ref) is ``g = \langle\cos\Theta\rangle = \alpha_1^1/3``
(MTL Eq. (4.92)), where ``\alpha_1^1`` is the first scattering-matrix expansion
coefficient.

### Scattering matrix

For a macroscopically isotropic and mirror-symmetric ensemble the scattering
matrix ``\mathbf{F}(\Theta)`` has the block-diagonal structure with up to ten
independent elements ``(\alpha_1\!-\!\alpha_4, \beta_1\!-\!\beta_6)``;
[`scattering_matrix`](@ref) evaluates it on a grid of scattering angles ``\Theta``.

## Methods for building ``T``

`TransitionMatrices.jl` provides several engines, all producing the same kind of
transition matrix:

| Method | Scatterers | Notes |
| --- | --- | --- |
| **Mie** (`bhmie`, `bhcoat`) | homogeneous & coated spheres | analytic; the reference for validation |
| **EBCM** (null-field) | axisymmetric (spheroid, cylinder, Chebyshev) | Waterman's method; fast, but the standard `Float64` integrands lose precision for elongated particles — use `stable` or extended precision |
| **IITM** (invariant imbedding) | axisymmetric *and* arbitrary (incl. prisms) | numerically robust where EBCM struggles; the only built-in route for non-axisymmetric shapes |
| **Sh-matrix** (moment separation) | axisymmetric | separates the geometry from ``(\lambda, m)`` so parameter **sweeps** reuse one quadrature |

Which one to pick, and when to enable `stable` or raise the precision, is covered
in [Choosing a solver](@ref). The literature behind each method is collected in
[Methods & references](@ref).
