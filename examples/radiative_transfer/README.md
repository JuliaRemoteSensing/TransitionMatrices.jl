# Radiative Transfer Workflow

This optional example adapts the WRF-driven dual-polarization radar prototype
contributed by [`@xiongyuup`](https://github.com/xiongyuup) in
[`PR #3`](https://github.com/JuliaRemoteSensing/TransitionMatrices.jl/pull/3).
It is a workflow example, not part of the `TransitionMatrices.jl` library API.

The workflow reads WRF-like NetCDF fields, precomputes a rain-drop scattering
table with `TransitionMatrices.jl`, integrates a gamma drop-size distribution at
each grid cell, and writes dual-pol radar fields to NetCDF.

## Setup

Run from the repository root:

```sh
julia --project=examples/radiative_transfer -e 'import Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

The dependencies live only in this example environment. They are not added to
the main package.

## Synthetic Smoke Test

The generator creates a tiny deterministic WRF-like file. It is only for
checking that I/O and the compute pipeline work; it is not validation data.

```sh
julia --project=examples/radiative_transfer \
  examples/radiative_transfer/generate_fake_wrf.jl \
  --output examples/radiative_transfer/fake_wrfout.nc
```

Then run the radar workflow:

```sh
JULIA_NUM_THREADS=4 julia --project=examples/radiative_transfer \
  examples/radiative_transfer/wrf_dualpol_radar.jl \
  --input examples/radiative_transfer/fake_wrfout.nc \
  --output examples/radiative_transfer/fake_dualpol.nc \
  --dsd double
```

The default diameter grid and IITM settings are intentionally small so the smoke
test finishes quickly. Use the CLI options to increase resolution for real
experiments.

## Inputs

The workflow expects these variables:

- `XLAT`, `XLONG`, `ZNU`, `XTIME`
- `QRAIN`, `QVAPOR`, `PB`, `P`, `T`
- `QNRAIN` for `--dsd double`

All 4-D microphysics and thermodynamic variables use
`west_east × south_north × bottom_top × Time` ordering.

## Outputs

The output NetCDF contains:

- `ZH`, `ZV`, `ZDR`, `RHOHV`, `KDP`
- copied `QRAIN` and derived/copied `QNRAIN`
- coordinate fields `XLAT`, `XLONG`, `ZNU`, `XTIME`

Dry or invalid cells are filled with `-999`.

## Notes

- `--dsd single` derives the number concentration from a fixed intercept
  parameter.
- `--dsd double` uses `QNRAIN` from the input file.
- A serialized scattering table cache is written by default to
  `scattering_table.jls`; delete it or change `--cache` when changing
  discretization settings beyond the checked metadata.
- Full-bin WRF microphysics is not implemented here. That remains a larger
  follow-up because it needs a different input contract and more careful unit
  handling.

## References

- H. Li, Y. Xiong, and Y. Chen, "Simulation of Complex Meteorological Target
  Echoes for Airborne Dual-Polarization Weather Radar Based on Invariant
  Imbedding T-Matrix," *IEEE Transactions on Geoscience and Remote Sensing*,
  62, 5105817, 2024.
- B. R. Brown, M. M. Bell, and A. J. Frambach, "Validation of Simulated
  Hurricane Drop Size Distributions Using Polarimetric Radar," *Geophysical
  Research Letters*, 42, 2016.
