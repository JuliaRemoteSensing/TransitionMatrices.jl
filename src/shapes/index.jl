export AbstractShape, AbstractHomogeneousShape, AbstractAxisymmetricShape, volume,
       volume_equivalent_radius, has_symmetric_plane

abstract type AbstractShape{T, CT} end
abstract type AbstractHomogeneousShape{T, CT} <: AbstractShape{T, CT} end
abstract type AbstractAxisymmetricShape{T, CT} <: AbstractHomogeneousShape{T, CT} end

volume(s::AbstractShape) = error("volume not implemented for $(typeof(s))")
volume_equivalent_radius(s::AbstractShape) = ∛(3.0 * volume(s) / (4.0 * π))
has_symmetric_plane(s::AbstractShape) = false

include("spheroid.jl")
include("chebyshev.jl")
