export TruncatedLinearProgrammingHall
export trunchall_dropdim
export reduce_same, reduce_sum, reduce_mean


"""
    TruncatedLinearProgrammingHall

A LinearProgrammingHall which drops one or more dimensions.

# Fields
* `component::Symbol`: The component defining a parameter or variable.
* `name::Symbol`: Either a parameter or variable name.
* `f::Vector{Float64}`: A vector of values, for every entry in the variable.
* `droppeddims::Vector{Symbol}`: A vector of dimensions to drop
"""
type TruncatedLinearProgrammingHall
    component::Symbol
    name::Symbol
    droppeddims::Vector{Symbol}
    f::Vector{Float64}
end

""" Return a TruncatedLinearProgrammingHall by dropping a dimension. """
function trunchall_dropdim(model::Model, hall::LinearProgrammingHall, dimension::Symbol, reduce::Function)
    dims = gettruncdims(model, hall.component, hall.name, [dimension])
    bigf = reshape(hall.f, getdims(model, hall.component, hall.name)...)
    f = Vector{Float64}(prod(dims))
    for ii in 1:length(f)
        index = untruncindex(ii, model, hall.component, hall.name, [dimension])
        f[ii] = reduce(bigf[index...])
    end

    TruncatedLinearProgrammingHall(hall.component, hall.name, [dimension], f)
end

"""
    TruncatedLinearProgrammingRoom

A LinearProgrammingRoom which drops one or more dimensions.

# Fields
* `varcomponent::Symbol`: The component defining the variable.
* `variable::Symbol`: The variable name.
* `vardropped::Vector{Symbol}`: A vector of dimensions to drop from the variable.
* `paramcomponent::Symbol`: The component defining the parameter.
* `parameter::Symbol`: The parameter name.
* `paramdropped::Vector{Symbol}`: A vector of dimensions to drop from the parameter.
* `A::SparseMatrixCSC{Float64, Int64}`: A sparse matrix of values, for every combination of the parameter and variable.
"""
type TruncatedLinearProgrammingRoom
    varcomponent::Symbol
    variable::Symbol
    vardropped::Vector{Symbol}
    paramcomponent::Symbol
    parameter::Symbol
    paramdropped::Vector{Symbol}
    A::SparseMatrixCSC{Float64, Int64}
end

function truncroom_droppardim(model::Model, room::LinearProgrammingRoom, dimension::Symbol, reduce::Function)
    vardims = getdims(model, room.varcomponent, room.variable)
    pardims = getdims(model, room.paramcomponent, room.parameter)
    truncpardims = gettruncdims(model, room.paramcomponent, room.parameter, [dimension])

    A = spzeros(prod(vardims), prod(truncpardims))

    for ii in 1:prod(vardims)
        bigf = reshape(room.A[ii, :], pardims...)
        for jj in 1:prod(truncpardims)
            index = untruncindex(jj, model, room.component, room.name, [dimension])
            value = reduce(bigf[index...])
            if value != 0
                A[ii, jj] = value
            end
        end
    end

    TruncatedLinearProgrammingRoom(room.varcomponent, room.variable, [],
                                   room.paramcomponent, room.parameter, [dimension], A)
end

"""
    TruncatedLinearProgrammingHouse

A LinearProgrammingHous which drops one or more dimnsions.

# Fields
* `model::Model`: The model containing all components
* `paramcomps::Vector{Symbol}`: The components defining each variable.
* `parameters::Vector{Symbol}`: The names of each variable.
* `paramdropped::Vector{Vector{Symbol}}`: By parameter, a vector of dimensions to drop.
* `constcomps::Vector{Symbol}`: The components defining each parameter.
* `constraints::Vector{Symbol}`: The names of each parameter.
* `constdropped::Vector{Vector{Symbol}}`: By constraint, a vector of dimensions to drop.
* `constdictionary::Dict{Symbol, Symbol}`: The names used in `constraints` must be unique, but can refer to the same parameter by including an entry in this disctionary mapping the unique name to the true variable name.
* `lowers::Vector{Float64}`: The lower bound for each parameter.
* `uppers::Vector{Float64}`: The upper bound for each parameter.
* `f::Vector{Float64}`: The derivative of the objective function for each parameter.
* `A::SparseMatrixCSC{Float64, Int64}`: Each row describes the derivatives of a given row for each parameter.
* `b::Vector{Float64}`: The maximum value for $A x$ for parameter values $x$.
"""
type LinearProgrammingHouse
    model::Model
    paramcomps::Vector{Symbol}
    parameters::Vector{Symbol}
    paramdropped::Vector{Vector{Symbol}}
    constcomps::Vector{Symbol}
    constraints::Vector{Symbol}
    constdropped::Vector{Vector{Symbol}}
    constdictionary::Dict{Symbol, Symbol}
    lowers::Vector{Float64}
    uppers::Vector{Float64}
    f::Vector{Float64} # length of parameters
    A::SparseMatrixCSC{Float64, Int64} # constraint rows, parameter columns
    b::Vector{Float64} # length of constraints
end

function trunchouse_droppardim(house::LinearProgrammingHouse, component::Symbol, parameter::Symbol, dimension::Symbol, reduce::Function)
    ## Start by getting varlengths
end

## Reduction functions

function reduce_same(arr)
    @assert all(arr .== arr[1])
    arr[1]
end

reduce_sum(arr) = sum(arr)
reduce_mean(arr) = mean(arr)

## Helpers

""" Return a vector of the dimension lengths, except for those dropped. """
function gettruncdims(model::Model, component::Symbol, name::Symbol, dropped::Vector{Symbol})
    dims = getdims(model, component, name)
    dimnames = getdimnames(model, component, name)

    convert(Vector{Int64}, [dims[ii] for ii in filter(ii -> !(dimnames[ii] in dropped), 1:length(dims))])
end

""" Return a vector of the dimension names, except for those droppd. """
function gettruncdimnames(model::Model, component::Symbol, name::Symbol, dropped::Vector{Symbol})
    dimnames = getdimnames(model, component, name)

    [dimname for dimname in filter(dimname -> !(dimname in dropped), dimnames)]
end

""" Like toindex for truncated elements, with colons in the dropped dimensions. """
function untruncindex(ii, model::Model, component::Symbol, name::Symbol, dropped::Vector{Symbol})
    dims = getdims(model, component, name)
    dimnames = getdimnames(model, component, name)

    # Dimensions, but with 1 for any truncated
    truncdims = convert(Vector{Int64}, [(dimnames[ii] in dropped ? 1 : dims[ii]) for ii in 1:length(dims)])
    truncindex = toindex(ii, truncdims)

    tuple([(dimnames[ii] in dropped ? (:) : truncindex[ii]) for ii in 1:length(dims)]...)
end
