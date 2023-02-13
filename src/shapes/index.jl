abstract type AbstractShape{T} end
abstract type AbstractHomogeneousShape{T} <: AbstractShape{T} end
abstract type AbstractAxisymmetricShape{T} <: AbstractHomogeneousShape{T} end
