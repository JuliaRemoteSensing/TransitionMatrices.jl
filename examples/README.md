# Examples

Runnable examples for `TransitionMatrices.jl`, as [Pluto.jl](https://plutojl.org)
notebooks. They use this `examples/` environment ([`Project.toml`](Project.toml)),
which includes `Pluto`, `Plots`, `ForwardDiff`, and develops the parent package.

## Running a notebook

```julia
julia> import Pkg; Pkg.activate("examples"); Pkg.instantiate()  # first time only
julia> using Pluto; Pluto.run()
```

Then open the notebook from the Pluto start page. Each notebook's first cell
activates this environment and develops `TransitionMatrices` from the repo root,
so it works from a fresh checkout (the first run installs dependencies).

## Notebooks

- [`shapes_gallery.jl`](shapes_gallery.jl) — **quick-start gallery**: define a
  spheroid, cylinder, Chebyshev particle, and N-fold prism; compute each
  T-matrix and its `Qsca`/`Qext`/`g`/`ω`; notes on which solver to use.
- [`solver_landscape.jl`](solver_landscape.jl) — **solver landscape &
  convergence**: the two-layer API (`EBCM`, `IITM`, `Iterative`, `stable`) on a
  high-aspect spheroid, showing classic EBCM diverging with `nₘₐₓ` while the
  stabilized path holds.
- [`angular_scattering.jl`](angular_scattering.jl) — **angular scattering &
  polarization**: the orientation-averaged scattering (Mueller) matrix vs angle
  (phase function `F₁₁`, linear polarization `-F₁₂/F₁₁`, `F₂₂/F₁₁`, `F₄₄/F₁₁`),
  plus the fixed-geometry amplitude/phase matrices.
- [`orientation_averaging.jl`](orientation_averaging.jl) — **orientation
  averaging**: Mishchenko's analytic `RandomOrientationTransitionMatrix` vs the
  numerical `orientation_average`, and the numerical average converging to it.
- [`spectral_sensitivity.jl`](spectral_sensitivity.jl) — **fast wavelength /
  refractive-index sensitivities** via the Sh-matrix moment-separation backend.
  `prepare_sh` runs the geometry quadrature once; differentiating the cheap
  reconstruction with `ForwardDiff` then gives exact `∂Cₛ𝚌ₐ/∂λ` and `∂Cₛ𝚌ₐ/∂mᵣ`
  across a whole spectrum at a large speedup over differentiating the
  from-scratch assembly. Includes plots of the spectra and their sensitivities,
  and a timing comparison.
- [`rain_radar.jl`](rain_radar.jl) — **rain radar observables**: adapt a
  radiative-transfer prototype contributed by
  [`@xiongyuup`](https://github.com/xiongyuup) in
  [PR #3](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/pull/3)
  into a self-contained example. It uses spheroidal IITM T-matrices, a
  rain-drop axis-ratio model, water/ice refractive-index fits, and a drop-size
  distribution to compute dual-pol radar moments and single- / multi-frequency
  brightness-temperature integrals without external WRF/NetCDF/MAT data.
