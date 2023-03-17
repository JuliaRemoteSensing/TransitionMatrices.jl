module TransitionMatrices

using Arblib
using Arblib: ArbLike, AcbLike, ArbVectorLike, AcbVectorLike, ArbMatrixLike, AcbMatrixLike
using ArbNumerics: ArbFloat, ArbReal, ArbComplex
using DoubleFloats: Double64
using FastGaussQuadrature: FastGaussQuadrature
using ForwardDiff: ForwardDiff
using GenericLinearAlgebra: cond, inv
using OffsetArrays: OffsetArray
using Quadmath: Quadmath, Float128, ComplexF128
using Rotations: Rotation, RotZYZ
using StaticArrays: SVector, SMatrix, SArray, @SVector, @SMatrix, @SArray
using TestItems: @testitem
using Wigxjpf: wig3jj, wig_table_init, wig_table_free, wig_temp_init, wig_temp_free

include("compat/index.jl")
include("special_functions/index.jl")
include("common/index.jl")
include("shapes/index.jl")

include("Mie/index.jl")
include("EBCM/index.jl")

# Various types of transition matrices
export AbstractTransitionMatrix, TransitionMatrix, RandomOrientationTransitionMatrix,
       AxisymmetricTransitionMatrix, MieTransitionMatrix

# Utility functions
export OrderDegreeIterator, rotate, amplitude_matrix, phase_matrix, scattering_matrix,
       orientation_average, scattering_cross_section, extinction_cross_section,
       absorption_cross_section,
       albedo, asymmetry_parameter,
       transition_matrix, transition_matrix_m, transition_matrix_m₀

# Utility functions (short names)
const calc_T = transition_matrix
const calc_S = amplitude_matrix
const calc_Z = phase_matrix
const calc_F = scattering_matrix
const calc_g = asymmetry_parameter
const calc_ω = albedo
const calc_Cext = extinction_cross_section
const calc_Csca = scattering_cross_section
const calc_Cabs = absorption_cross_section

export calc_T, calc_S, calc_Z, calc_F, calc_ω,
       calc_Csca, calc_Cext, calc_Cabs

# Shape related exports
export AbstractShape, AbstractHomogeneousShape, AbstractAxisymmetricShape, volume,
       volume_equivalent_radius, has_symmetric_plane, Spheroid, Cylinder, Chebyshev

# Re-exports
export RotZYZ, Double64, Float128, ComplexF128, Arb, Acb, ArbFloat, ArbReal, ArbComplex

end
