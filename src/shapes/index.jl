abstract type AbstractShape{T, CT} end
abstract type AbstractAxisymmetricShape{T, CT} <: AbstractShape{T, CT} end
abstract type AbstractNFoldShape{N, T, CT} <: AbstractShape{T, CT} end

volume(s::AbstractShape) = 4 // 3 * π * volume_equivalent_radius(s)^3

# The following functions should be implemented for each concrete shape type
volume_equivalent_radius(::AbstractShape{T}) where {T} = zero(T)
refractive_index(::AbstractShape{T, CT}, x) where {T, CT} = one(CT)
has_symmetric_plane(::AbstractShape) = false

# Axisymmetric shapes

include("spheroid.jl")
include("cylinder.jl")
include("chebyshev.jl")

# NFold shapes

include("prism.jl")
