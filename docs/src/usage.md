```@meta
CurrentModule = TransitionMatrices
```

# Usage

## Define a shape

### Spheroid

The following defines a spheroid with semi-major axis `a=1.0`, semi-minor axis `c=0.5`, and a refractive index of `m=1.311`.

```julia
spheroid = Spheroid{Float64, ComplexF64}(1.0, 0.5, 1.311);
```

The type specification `{Float64, ComplexF64}` can be omitted, and in that case, the real and complex types used in the calculations will be inferred from the input parameters.

!!! note

    In the example above, if `{Float64, ComplexF64}` is ommited, do remember to use `complex(1.311)` instead of 1.311 to denote that the refractive index is a complex value. Otherwise errors would occur in later calculations.

### Cylinder

The following defines a cylinder with radius `r=1.0`, height `h=2.0`, and a refractive index of `m=1.5+0.01im`. Note that the type specification `{Float128, ComplexF128}` is used to specify the real and complex types used in the calculations.

```julia
cylinder = Cylinder{Float128, ComplexF128}(1.0, 2.0, 1.5+0.01im);
```

### Chebyshev particle

The following defines a Chebyshev particle with radius `r₀=1.0`, deformation coefficient `ε=0.1`, order `n=4`, and a refractive index of `m=1.5+0.01im`. Note that the type specification `{Arb, Acb}` is used to specify the real and complex types used in the calculations.

```julia
chebyshev = Chebyshev{Arb, Acb}(1.0, 0.1, 4, 1.5+0.01im);
```

## Calculate the T-Matrix

After defining a shape, the T-Matrix can be calculated using the `transition_matrix` function. It has an alias `calc_T` for convenience.

Take the spheroid defined above as an example, and suppose the wavelength is `λ=2π`. The T-Matrix can be calculated iteratively as follows:

```julia
𝐓 = calc_T(spheroid, 2π)
```

To see the process, you can turn on debugging by setting `ENV["JULIA_DEBUG"] = "TransitionMatrices"`, and then run the above code.

There are many optional keyword arguments, and you can refer to the detailed documentation of the [`transition_matrix`](@ref).

You can also calculate the T-Matrix directly by specifying the truncation order and number of Gauss-Legendre quadrature points:

```julia
𝐓 = transition_matrix(spheroid, 2π, 10, 100)
```

## Post-processing

After getting the T-Matrix, you can calculate the far-field scattering properties using the following functions:

- [`amplitude_matrix`](@ref), a.k.a. `calc_S`
- [`phase_matrix`](@ref), a.k.a. `calc_Z`
- [`extinction_cross_section`](@ref), a.k.a. `calc_Cext`
- [`scattering_cross_section`](@ref), a.k.a. `calc_Csca`
- [`absorption_cross_section`](@ref), a.k.a. `calc_Cabs`
- [`albedo`](@ref), a.k.a. `calc_ω`
- [`asymmetry_parameter`](@ref), a.k.a. `calc_g`

And the orientation-averaged scattering matrix:

- [`scattering_matrix`](@ref), a.k.a. `calc_F`

## Automatic differentiation

`TransitionMatrices.jl` supports automatic differentiation using [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl).

Here is an example of calculating the gradient of the scattering cross section with respect to `a`, `c`, `mᵣ`, `mᵢ` and `λ`:

```jldoctest
using ForwardDiff, TransitionMatrices

function f(x)
    s = Spheroid(x[1], x[2], complex(x[3], x[4]))
    T₀ = TransitionMatrices.transition_matrix(s, x[5], 10, 100)
    Csca = calc_Csca(T₀)
end

gradient = ForwardDiff.gradient(f, [2.0, 3.0, 1.311, 0.02, 2π])

# output

5-element Vector{Float64}:
  19.872009666582716
   5.740718904091895
  89.1305069099317
 -45.845687597163966
  -9.06644850667507
```
