"""
    AbstractLinearizationBackend

Abstract supertype for analytical or linearized differentiation backends.
Concrete backends declare which solvers can provide a Jacobian instead of
falling back to numerical differentiation.
"""
abstract type AbstractLinearizationBackend end

function _linearization_property(source, name::Symbol; default = nothing)
    isnothing(source) && return default
    return hasproperty(source, name) ? getproperty(source, name) : default
end

function _linearization_variables_supported(vars, canonical)
    all(var -> var in canonical, vars) || return false
    return length(Set(vars)) == length(vars)
end

"""
    MieLinearization()

Backend marker for future analytical derivatives of Mie transition matrices.
"""
struct MieLinearization <: AbstractLinearizationBackend end

"""
    EBCMLinearization()

Backend marker for linearized EBCM transition matrices. The current supported
slice uses analytical `P`/`U` derivatives for fixed-order, fixed-quadrature
Spheroid, Chebyshev, and Cylinder problems, then propagates those derivatives
through the analytical EBCM matrix identity.
"""
struct EBCMLinearization <: AbstractLinearizationBackend end

"""
    IITMLinearization([variant])

Backend marker for future analytical derivatives of IITM transition matrices.
`variant` distinguishes `:axisymmetric`, `:nfold`, `:arbitrary`, or `:auto`.
"""
struct IITMLinearization <: AbstractLinearizationBackend
    variant::Symbol
end

IITMLinearization() = IITMLinearization(:auto)

"""
    LinearizationProblem(rebuild, x; variables)
    LinearizationProblem(x; variables, rebuild)

Describe a real parameter vector and the function that rebuilds the physical
scattering input from that vector. `variables` must name each entry of `x`.

This mirrors the existing numerical differentiation pattern where users write
`x -> shape, wavelength, config`, while giving analytical backends a stable
problem object to dispatch on.
"""
struct LinearizationProblem{TX <: AbstractVector, F, V}
    x::TX
    variables::V
    rebuild::F
end

function LinearizationProblem(rebuild::F, x::AbstractVector; variables) where {F}
    vars = Tuple(variables)
    length(vars) == length(x) ||
        throw(ArgumentError("linearization variable count must match parameter count"))
    all(v -> v isa Symbol, vars) ||
        throw(ArgumentError("linearization variables must be symbols"))

    return LinearizationProblem{typeof(x), F, typeof(vars)}(x, vars, rebuild)
end

function LinearizationProblem(x::AbstractVector; variables, rebuild)
    return LinearizationProblem(rebuild, x; variables)
end

variables(problem::LinearizationProblem) = problem.variables
Base.length(problem::LinearizationProblem) = length(problem.x)
Base.getindex(problem::LinearizationProblem, i) = problem.x[i]
rebuild(problem::LinearizationProblem, x = problem.x) = problem.rebuild(x)

"""
    LinearizationResult(value, jacobian, variables; metadata=(;))

Container for a computed value and its Jacobian with respect to a
`LinearizationProblem` parameter vector.
"""
struct LinearizationResult{V, J, Vars, M}
    value::V
    jacobian::J
    variables::Vars
    metadata::M
end

function LinearizationResult(value, jacobian, variables; metadata = (;))
    vars = Tuple(variables)
    return LinearizationResult{typeof(value), typeof(jacobian), typeof(vars),
        typeof(metadata)}(value, jacobian, vars, metadata)
end

variables(result::LinearizationResult) = result.variables

"""
    derivative(result, i)
    derivative(result, variable)

Return the derivative slice for the `i`th variable, or for the named
`variable`. If the Jacobian stores parameter derivatives along its last
dimension, the returned object is a view of that slice.
"""
function _linearization_variable_index(vars, variable::Symbol)
    i = findfirst(isequal(variable), vars)
    isnothing(i) && throw(ArgumentError("unknown linearization variable: $variable"))
    return i
end

function derivative(result::LinearizationResult, i::Integer)
    1 <= i <= length(result.variables) || throw(BoundsError(result.variables, i))
    J = result.jacobian

    if result.value isa AbstractArray && J isa AbstractArray &&
       ndims(J) == ndims(result.value) + 1 && size(J, ndims(J)) == length(result.variables)
        return selectdim(J, ndims(J), i)
    end

    return J[i]
end

function derivative(result::LinearizationResult, variable::Symbol)
    return derivative(result, _linearization_variable_index(result.variables, variable))
end

"""
    LinearizationSupport(supported, reason)

Capability record returned by [`supports_linearization`](@ref).
"""
struct LinearizationSupport
    supported::Bool
    reason::String
end

Base.Bool(support::LinearizationSupport) = support.supported

"""
    UnsupportedLinearization

Exception thrown when the unified linearization API is asked for an analytical
Jacobian that has not been implemented for the chosen problem/backend/output.
"""
struct UnsupportedLinearization <: Exception
    backend::AbstractLinearizationBackend
    output::Symbol
    reason::String
end

function Base.showerror(io::IO, err::UnsupportedLinearization)
    print(io, "Unsupported linearization for output :", err.output,
        " with backend ", typeof(err.backend), ": ", err.reason)
end

"""
    supports_linearization(problem, backend; output=:transition_matrix, config=nothing)

Return a [`LinearizationSupport`](@ref) record describing whether an analytical
linearization implementation is available for the requested output.
"""
function supports_linearization(::LinearizationProblem,
        backend::AbstractLinearizationBackend;
        output::Symbol = :transition_matrix,
        config = nothing)
    return LinearizationSupport(false,
        "no analytical linearization implementation is registered")
end

"""
    linearize_transition_matrix(problem, backend; config=nothing)

Compute a transition matrix and its Jacobian with respect to `problem`.
Backends that do not implement analytical derivatives throw
[`UnsupportedLinearization`](@ref).
"""
function linearize_transition_matrix(problem::LinearizationProblem,
        backend::AbstractLinearizationBackend;
        config = nothing)
    support = supports_linearization(problem, backend; output = :transition_matrix, config)
    throw(UnsupportedLinearization(backend, :transition_matrix, support.reason))
end

"""
    linearize_observable(observable, problem, backend; config=nothing)

Compute a scattering observable and its Jacobian with respect to `problem`.
The default method only defines the unsupported boundary; analytical backends
should add methods that reuse `linearize_transition_matrix`.
"""
function linearize_observable(observable, problem::LinearizationProblem,
        backend::AbstractLinearizationBackend; config = nothing)
    output = observable isa Symbol ? observable : :observable
    support = supports_linearization(problem, backend; output, config)
    throw(UnsupportedLinearization(backend, output, support.reason))
end
