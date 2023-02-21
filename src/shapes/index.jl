abstract type AbstractShape{T, CT} end
abstract type AbstractHomogeneousShape{T, CT} <: AbstractShape{T, CT} end
abstract type AbstractAxisymmetricShape{T, CT} <: AbstractHomogeneousShape{T, CT} end

volume(s::AbstractShape) = 4 // 3 * π * volume_equivalent_radius(s)^3

# The following function should be implemented for each concrete shape type
# volume_equivalent_radius(s::AbstractShape)

has_symmetric_plane(s::AbstractShape) = false

include("spheroid.jl")
include("cylinder.jl")
include("chebyshev.jl")
