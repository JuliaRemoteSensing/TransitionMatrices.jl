module TransitionMatrices

using Arblib
using Arblib: ArbLike, AcbLike, ArbVectorLike, AcbVectorLike, ArbMatrixLike, AcbMatrixLike
using DoubleFloats: Double64
using FastGaussQuadrature: FastGaussQuadrature
using LinearAlgebra: cond
using OffsetArrays: OffsetArray
using Quadmath: Float128
using Rotations: Rotation, RotZYZ
using StaticArrays: SVector, SMatrix, SArray, @SVector, @SMatrix, @SArray
using TestItems: @testitem

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
export OrderDegreeIterator, rotate, amplitude_matrix, phase_matrix,
       orientation_average, transition_matrix,
       transition_matrix_m, transition_matrix_m₀

# Shape related exports
export AbstractShape, AbstractHomogeneousShape, AbstractAxisymmetricShape, volume,
       volume_equivalent_radius, has_symmetric_plane, Spheroid, Chebyshev

# Re-exports
export RotZYZ, Double64, Float128, Arb, Acb

end
