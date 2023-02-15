abstract type AbstractShape{T, CT} end
abstract type AbstractHomogeneousShape{T, CT} <: AbstractShape{T, CT} end
abstract type AbstractAxisymmetricShape{T, CT} <: AbstractHomogeneousShape{T, CT} end

include("spheroid.jl")
