using MathProgBase
using DataFrames

import Mimi.metainfo
import Base.*, Base.-, Base.+, Base.max

export LinearProgrammingHall, LinearProgrammingShaft, LinearProgrammingRoom, LinearProgrammingHouse
export hallsingle, hall_relabel, hallvalues
export shaftsingle, shaftvalues
export roomdiagonal, roomsingle, roomempty, roomintersect, room_relabel, room_relabel_parameter
export varsum, fromindex
export setobjective!, setconstraint!, setconstraintoffset!, getconstraintoffset, clearconstraint!, setlower!, setupper!, addparameter!, addconstraint!
## Debugging
export getroom, getnonzerorooms
export gethouse, constraining, houseoptimize, summarizeparameters, findinfeasiblepair, varlengths, getconstraintsolution, getparametersolution

include("metamimi.jl")

# A hallway is a vector of variables
"""
    LinearProgrammingHall

A vector of values, used either to describes the objective function or
the constraint offsets.

# Fields
* `component::Symbol`: The component defining a parameter or variable.
* `name::Symbol`: Either a parameter or variable name.
* `f::Vector{Float64}`: A vector of values, for every entry in the variable.
"""
type LinearProgrammingHall
    component::Symbol
    name::Symbol
    f::Vector{Float64}
end

"Construct a hall, calling gen with each index."
function hallsingle(model::Model, component::Symbol, name::Symbol, gen::Function)
    LinearProgrammingHall(component, name, vectorsingle(getdims(model, component, name), gen))
end

function hallvalues(model::Model, component::Symbol, name::Symbol, values::Array{Float64})
    gen(inds...) = values[inds...]
    hallsingle(model, component, name, gen)
end

"""Connect a derivative to another component: change the variable component and name to another component."""
function hall_relabel(hall::LinearProgrammingHall, from::Symbol, tocomponent::Symbol, toname::Symbol)
    @assert hall.name == from "Name mismatch in hall_relabel: $(hall.name) <> $from"

    LinearProgrammingHall(tocomponent, toname, hall.f)
end

function -(hall::LinearProgrammingHall)
    LinearProgrammingHall(hall.component, hall.name, -hall.f)
end

"""
    hall1 + hall2

Construct a new LinearProgrammingHall which is the element-wise
sum of the values in hall1 and hall2, which must refer to
the same variable.
"""
function +(hall1::LinearProgrammingHall, hall2::LinearProgrammingHall; skipnamecheck=false)
    if !skipnamecheck
        @assert hall1.name == hall2.name "Hall + Hall name mismatch: $(hall1.name) <> $(hall2.name); use hall_relabel?"
    end
    LinearProgrammingHall(hall1.component, hall1.name, hall1.f + hall2.f)
end

"""
    hall1 - hall2

Construct a new LinearProgrammingHall which is the element-wise
difference in the values in hall1 and hall2, which must refer to
the same variable.
"""
function -(hall1::LinearProgrammingHall, hall2::LinearProgrammingHall; skipnamecheck=false)
    if !skipnamecheck
        @assert hall1.name == hall2.name "Hall - Hall name mismatch: $(hall1.name) <> $(hall2.name); use hall_relabel?"
    end
    LinearProgrammingHall(hall1.component, hall1.name, hall1.f - hall2.f)
end

"""
    max(hall1, hall2)

Construct a new LinearProgrammingHall which is the element-wise
maximum value of each entry in hall1 and hall2, which must refer to
the same variable.
"""
function max(hall1::LinearProgrammingHall, hall2::LinearProgrammingHall; skipnamecheck=false)
    if !skipnamecheck
        @assert hall1.name == hall2.name "Hall - Hall name mismatch: $(hall1.name) <> $(hall2.name); use hall_relabel?"
    end
    LinearProgrammingHall(hall1.component, hall1.name, max(hall1.f, hall2.f))
end

# A shaft is a transpose of a hallway, used for parameters
"""
    LinearProgrammingShaft

A vector of values, used either to describes the objective function or
the constraint offsets.  This acts as the transpose of a Hall,
although this distinction is only important for multiplying a Room by
Shaft, which then returns a Hall.

# Fields
* `component::Symbol`: The component defining a parameter or variable.
* `name::Symbol`: Either a parameter or variable name.
* `x::Vector{Float64}`: A vector of values, for every entry in the variable.
"""
type LinearProgrammingShaft
    component::Symbol
    name::Symbol
    x::Vector{Float64}
end

function shaftsingle(model::Model, component::Symbol, name::Symbol, gen::Function)
    LinearProgrammingShaft(component, name, vectorsingle(getdims(model, component, name), gen))
end

function shaftvalues(model::Model, component::Symbol, name::Symbol, values::Array{Float64})
    gen(inds...) = values[inds...]
    shaftsingle(model, component, name, gen)
end

function -(shaft::LinearProgrammingShaft)
    LinearProgrammingShaft(shaft.component, shaft.name, -shaft.x)
end

function +(shaft1::LinearProgrammingShaft, shaft2::LinearProgrammingShaft; skipnamecheck=false)
    if !skipnamecheck
        @assert shaft1.name == shaft2.name "Shaft + Shaft name mismatch: $(shaft1.name) <> $(shaft2.name); use shaft_relabel?"
    end
    LinearProgrammingShaft(shaft1.component, shaft1.name, shaft1.x + shaft2.x)
end

"""
    LinearProgrammingRoom

A matrix of values, used to describe a part of the linear programming
problem constraint matrix (a gradient matrix).  The rows of the matrix
correspond to variable indices, and the columns correspond to
parameter indices.  The values are stored as a sparse matrix, since
most parameter indices are assumed to not interact with most variable
indices.

# Fields
* `varcomponent::Symbol`: The component defining the variable.
* `variable::Symbol`: The variable name.
* `paramcomponent::Symbol`: The component defining the parameter.
* `parameter::Symbol`: The parameter name.
* `A::SparseMatrixCSC{Float64, Int64}`: A sparse matrix of values, for every combination of the parameter and variable.
"""
type LinearProgrammingRoom
    varcomponent::Symbol
    variable::Symbol
    paramcomponent::Symbol
    parameter::Symbol
    A::SparseMatrixCSC{Float64, Int64}
end

doc"""
    roomdiagonal(model, component, variable, parameter, gen)

Fill in just the diagonal, assuming $\frac{dv_i}{dp_j} = 0$, if $i <>
j$.  Requires that the dimensions of `variable` and `parameter` are
the same.  `gen` called for each combination of indices.

$$\begin{array}{ccc} p_1 & p_2 & \cdots \end{array}$$
$$\begin{array}{c}
    v_1 \\
    v_2 \\
    \vdots
    \end{array}\left(\begin{array}{ccc}
        g(1) & 0 & \cdots \\
        0 & g(2) & \cdots \\
        \vdots & \vdots & \ddots
        \end{array}\right)$$

# Arguments
* `model::Model`: The model containing `component`.
* `component::Symbol`: The component containing `variable` and `parameter`.
* `variable::Symbol`: The outcome variable, corresponding to matrix rows.
* `parameter::Symbol`: The optimization parameter, corresponding to matrix columns.
* `gen::Function`: The function generating gradient values.
"""
function roomdiagonal(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, gen::Function)
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)
    @assert dimsvar == dimspar "Variable and parameter in roomdiagonal do not have the same dimensions: $dimsvar <> $dimspar"
    LinearProgrammingRoom(component, variable, component, parameter, matrixdiagonal(dimsvar, gen))
end

"""
Fill in every element
"""
function roomsingle(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, gen::Function)
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)
    LinearProgrammingRoom(component, variable, component, parameter, matrixsingle(dimsvar, dimspar, gen))
end

"""
    roomintersect(model, component, variable1, variable2, gen)

Generate a room at the intersection of two variables.  Variable1 is
the constraint, and variable2 has the optimization parameter.

Call `gen` for each shared index, passing an array to be filled for unshared indexes.
This version assumes both variables come from the same component.

See `matrixintersect` for the matrix generation logic.
"""
function roomintersect(model::Model, component::Symbol, variable1::Symbol, variable2::Symbol, gen::Function)
    dims1 = getdims(model, component, variable1)
    dims2 = getdims(model, component, variable2)
    dimnames1 = getdimnames(model, component, variable1)
    dimnames2 = getdimnames(model, component, variable2)
    LinearProgrammingRoom(component, variable1, component, variable2, matrixintersect(dims1, dims2, dimnames1, dimnames2, gen))
end

"""
    roomintersect(model, component, variable1, component2, variable2, gen)

Generate a room at the intersection of two variables.

Call `gen` for each shared index, passing an array to be filled for unshared indexes.
This version allows the variables to come from different component.

See `matrixintersect` for the matrix generation logic.
"""
function roomintersect(model::Model, component1::Symbol, variable1::Symbol, component2::Symbol, variable2::Symbol, gen::Function)
    dims1 = getdims(model, component1, variable1)
    dims2 = getdims(model, component2, variable2)
    dimnames1 = getdimnames(model, component1, variable1)
    dimnames2 = getdimnames(model, component2, variable2)
    LinearProgrammingRoom(component1, variable1, component2, variable2, matrixintersect(dims1, dims2, dimnames1, dimnames2, gen))
end

"""
Create a room to be filled in later
"""
function roomempty(model::Model, component::Symbol, variable::Symbol, parameter::Symbol)
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)
    LinearProgrammingRoom(component, variable, component, parameter, matrixempty(dimsvar, dimspar))
end

"""Connect a gradient to another component: change the variable component and name to another component."""
function room_relabel(room::LinearProgrammingRoom, from::Symbol, tocomponent::Symbol, toname::Symbol)
    @assert room.variable == from "Variable name mismatch in room_relabel: $(room.variable) <> $from"

    LinearProgrammingRoom(tocomponent, toname, room.paramcomponent, room.parameter, room.A)
end

"""Connect a gradient to another component: change the variable component and name to another component."""
function room_relabel_parameter(room::LinearProgrammingRoom, from::Symbol, tocomponent::Symbol, toname::Symbol)
    @assert room.parameter == from "Parameter name mismatch in room_relabel_parameter: $(room.variable) <> $from"

    LinearProgrammingRoom(room.varcomponent, room.variable, tocomponent, toname, room.A)
end

"""
    Construct a new LinearProgrammingHall, with each entry as a weighted sum of the partials for a given variable.
        This is equivalent to translating a constraint by a functional relationship.
"""
function *(hall::LinearProgrammingHall, room::LinearProgrammingRoom; skipnamecheck=false)
    if !skipnamecheck
        @assert hall.name == room.variable "Hall * Room name mismatch: $(hall.name) <> $(room.variable); use room_relabel?"
    end
    LinearProgrammingHall(room.paramcomponent, room.parameter, vec(transpose(hall.f) * room.A))
end

# Returns a hall, since it then contains variable/constraint information
function *(room::LinearProgrammingRoom, shaft::LinearProgrammingShaft; skipnamecheck=false)
    if !skipnamecheck
        @assert shaft.name == room.parameter "Room * Shaft name mismatch: $(room.parameter) <> $(hall.name); use room_relabel?"
    end
    LinearProgrammingHall(room.varcomponent, room.variable, room.A * shaft.x)
end


function *(room1::LinearProgrammingRoom, room2::LinearProgrammingRoom; skipnamecheck=false)
    if !skipnamecheck
        @assert room1.parameter == room2.variable "Room * Room name mismatch: $(room1.parameter) <> $(room2.variable); use room_relabel?"
    end
    LinearProgrammingRoom(room1.varcomponent, room1.variable, room2.paramcomponent, room2.parameter, room1.A * room2.A)
end

function +(room1::LinearProgrammingRoom, room2::LinearProgrammingRoom; skipnamecheck=false)
    if !skipnamecheck
        @assert room1.parameter == room2.parameter "Room + Room parameter name mismatch: $(room1.parameter) <> $(room2.parameter); use room_relabel?"
        @assert room1.variable == room2.variable "Room + Room variable name mismatch: $(room1.variable) <> $(room2.variable); use room_relabel?"
    end
    if size(room1.A)[1] == 0 || size(room1.A)[2] == 0
        LinearProgrammingRoom(room1.varcomponent, room1.variable, room1.paramcomponent, room1.parameter, spzeros(size(room1.A)...))
    else
        LinearProgrammingRoom(room1.varcomponent, room1.variable, room1.paramcomponent, room1.parameter, room1.A + room2.A)
    end
end

function -(room1::LinearProgrammingRoom, room2::LinearProgrammingRoom; skipnamecheck=false)
    if !skipnamecheck
        @assert room1.parameter == room2.parameter "Room - Room parameter name mismatch: $(room1.parameter) <> $(room2.parameter); use room_relabel?"
        @assert room1.variable == room2.variable "Room - Room variable name mismatch: $(room1.variable) <> $(room2.variable); use room_relabel?"
    end
    if size(room1.A)[1] == 0 || size(room1.A)[2] == 0
        LinearProgrammingRoom(room1.varcomponent, room1.variable, room1.paramcomponent, room1.parameter, spzeros(size(room1.A)...))
    else
        LinearProgrammingRoom(room1.varcomponent, room1.variable, room1.paramcomponent, room1.parameter, room1.A - room2.A)
    end
end

function -(room::LinearProgrammingRoom)
    LinearProgrammingRoom(room.varcomponent, room.variable, room.paramcomponent, room.parameter, -room.A)
end

"""
Construct a hall from a room by summing over rows.
        This is equivalent to a constraint for the sum of the variables in the original room.
"""
function varsum(room::LinearProgrammingRoom)
    LinearProgrammingHall(room.paramcomponent, room.parameter, vec(sum(room.A, 1)))
end

doc"""
    LinearProgrammingHouse

The full description of a linear programming problem, including all
its variables, parameters, the constraints, and the objective.

The linear programming that is solved is always:
$$\max f' x$$
$$A x \le b$$
$$x_{lower} \le x \le x_{upper}$$

# Fields
* `model::Model`: The model containing all components
* `paramcomps::Vector{Symbol}`: The components defining each variable.
* `parameters::Vector{Symbol}`: The names of each variable.
* `constcomps::Vector{Symbol}`: The components defining each parameter.
* `constraints::Vector{Symbol}`: The names of each parameter.
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
    constcomps::Vector{Symbol}
    constraints::Vector{Symbol}
    constdictionary::Dict{Symbol, Symbol}
    lowers::Vector{Float64}
    uppers::Vector{Float64}
    f::Vector{Float64} # length of parameters
    A::SparseMatrixCSC{Float64, Int64} # constraint rows, parameter columns
    b::Vector{Float64} # length of constraints
end

function LinearProgrammingHouse(model::Model, paramcomps::Vector{Symbol}, parameters::Vector{Symbol}, constcomps::Vector{Symbol}, constraints::Vector{Symbol}, constdictionary::Dict{Symbol, Symbol}=Dict{Symbol, Symbol}())
    paramlen = sum(varlengths(model, paramcomps, parameters))
    variablelen = sum(varlengths(model, constcomps, constraints, constdictionary))
    A = spzeros(variablelen, paramlen)
    f = zeros(paramlen)
    b = zeros(variablelen)
    LinearProgrammingHouse(model, paramcomps, parameters, constcomps, constraints, constdictionary, zeros(length(f)), Inf * ones(length(f)), f, A, b)
end

function addparameter!(house::LinearProgrammingHouse, component::Symbol, parameter::Symbol)
    paramlen = prod(getdims(house.model, component, parameter))

    append!(house.paramcomps, [component])
    append!(house.parameters, [parameter])
    append!(house.lowers, zeros(paramlen))
    append!(house.uppers, Inf * ones(paramlen))
    append!(house.f, zeros(paramlen))
    house.A = [house.A spzeros(length(house.b), paramlen)]
end

function addconstraint!(house::LinearProgrammingHouse, component::Symbol, constraint::Symbol, variable::Symbol)
    varlen = prod(getdims(house.model, component, variable))

    append!(house.constcomps, [component])
    append!(house.constraints, [constraint])
    house.constdictionary[constraint] = variable
    append!(house.b, zeros(varlen))
    house.A = [house.A; spzeros(varlen, length(house.f))]
end

function setobjective!(house::LinearProgrammingHouse, hall::LinearProgrammingHall)
    @assert hall.name in house.parameters "$(hall.name) not a known parameter"
    kk = findfirst((house.paramcomps .== hall.component) & (house.parameters .== hall.name))
    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    @assert length(hall.f) == paramspans[kk] "Length of parameter $(hall.name) unexpected: $(length(hall.f)) <> $(paramspans[kk])"
    before = sum(paramspans[1:kk-1])
    for ii in 1:paramspans[kk]
        #@assert house.f[before+ii] == 0 "Overwrite existing gradient in setobjective!"
        house.f[before+ii] = hall.f[ii]
    end
end

function setconstraint!(house::LinearProgrammingHouse, room::LinearProgrammingRoom)
    kk = findfirst((house.paramcomps .== room.paramcomponent) & (house.parameters .== room.parameter))
    @assert kk > 0 "$(room.paramcomponent).$(room.parameter) not a known parameter"
    ll = findfirst((house.constcomps .== room.varcomponent) & (house.constraints .== room.variable))
    @assert ll > 0 "$(room.varcomponent).$(room.variable) not a known variable"

    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)
    @assert size(room.A, 1) == constspans[ll] "Length of variable $(room.variable) unexpected: $(size(room.A, 1)) <> $(constspans[ll])"
    @assert size(room.A, 2) == paramspans[kk] "Length of parameter $(room.parameter) unexpected: $(size(room.A, 2)) <> $(paramspans[kk])"

    if size(room.A, 2) == 0 # expressions below mishandle this case
        return
    end

    parambefore = sum(paramspans[1:kk-1])
    constbefore = sum(constspans[1:ll-1])

    house.A[constbefore+1:constbefore+constspans[ll], parambefore+1:parambefore+paramspans[kk]] = room.A
end

function getroom(house::LinearProgrammingHouse, varcomponent::Symbol, variable::Symbol, paramcomponent::Symbol, parameter::Symbol)
    @assert parameter in house.parameters "$parameter not a known parameter"
    @assert variable in house.constraints "$variable not a known variable"

    # Determine where to index into vector
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)
    rowkk = findfirst((house.constcomps .== varcomponent) & (house.constraints .== variable))
    rowbefore = sum(constspans[1:rowkk-1])

    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    colkk = findfirst((house.paramcomps .== paramcomponent) & (house.parameters .== parameter))
    colbefore = sum(paramspans[1:colkk-1])

    # Collect the offset
    subA = house.A[(rowbefore + 1):(rowbefore + constspans[rowkk]), (colbefore + 1):(colbefore + paramspans[colkk])]

    LinearProgrammingRoom(varcomponent, variable, paramcomponent, parameter, subA)
end

"""Return a list of all rooms (parameter+component and variable+component) that are related by the derivative matrix.
"""
function getnonzerorooms(house::LinearProgrammingHouse)
    df = DataFrame(varcomponent=Symbol[], variable=Symbol[], paramcomponent=Symbol[], parameter=Symbol[])
    for kk in 1:length(house.paramcomps)
        for ll in 1:length(house.constcomps)
            room = getroom(house, house.constcomps[ll], house.constraints[ll], house.paramcomps[kk], house.parameters[kk])
            if !isempty(nonzeros(room.A))
                push!(df, [house.constcomps[ll], house.constraints[ll], house.paramcomps[kk], house.parameters[kk]])
            end
        end
    end

    df
end


"""
Set offset values from a `LinearProgrammingHall`.  See the other `setconstraintoffset!`
"""
function setconstraintoffset!(house::LinearProgrammingHouse, hall::LinearProgrammingHall)
    setconstraintoffset!(house, hall.component, hall.name, hall.f)
end

"""
    setconstraintoffset!(house, component, variable, f)

Set offset values within a `LinearProgrammingHouse`.

# Arguments
* `house::LinearProgrammingHouse`: The house to set values within.
* `component::Symbol`: The component for the constraint variable.
* `variable::Symbol`: The variable for the constraint.
* `f::Vector{Float64}`: The values, with all dimensions collapsed into a single vector.
"""
function setconstraintoffset!(house::LinearProgrammingHouse, component::Symbol, variable::Symbol, f::Vector{Float64})
    @assert variable in house.constraints "$(variable) not a known variable"
    kk = findfirst((house.constcomps .== component) & (house.constraints .== variable))
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)
    @assert length(f) == constspans[kk] "Length of parameter $(variable) unexpected: $(length(f)) <> $(constspans[kk])"
    before = sum(constspans[1:kk-1])
    for ii in 1:constspans[kk]
        #@assert house.b[before+ii] == 0 "Overwrite existing gradient in setobjective!"
        house.b[before+ii] = f[ii]
    end
end

"""
    getconstraintoffset(house, component, variable, reshp)

Return the values for a constraint, optionally reshaped to the original dimensions.

# Arguments
* `house::LinearProgrammingHouse`: The house from which to get the values.
* `component::Symbol`: The component for the constraint variable.
* `variable::Symbol`: The variable for the constraint.
* `reshp::Bool`: Should it be reshaped to the original variable dimensions? (default: false)
"""
function getconstraintoffset(house::LinearProgrammingHouse, component::Symbol, variable::Symbol; reshp::Bool=false)
    @assert variable in house.constraints "$variable not a known variable"

    # Determine where to index into vector
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)
    kk = findfirst((house.constcomps .== component) & (house.constraints .== variable))
    before = sum(constspans[1:kk-1])

    # Collect the offset
    constraintoffset = house.b[(before + 1):(before + constspans[kk])]

    if reshp
        # Reshape to original dimensions
        reshape(constraintoffset, getdims(house.model, component, variable)...)
    else
        constraintoffset
    end
end

function clearconstraint!(house::LinearProgrammingHouse, component::Symbol, variable::Symbol)
    @assert variable in house.constraints "$(variable) not a known variable"

    ll = findfirst((house.constcomps .== component) & (house.constraints .== variable))
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)

    constbefore = sum(constspans[1:ll-1])

    house.A[constbefore+1:constbefore+constspans[ll], :] = 0
    house.b[constbefore+1:constbefore+constspans[ll]] = 0
end

function setlower!(house::LinearProgrammingHouse, hall::LinearProgrammingHall)
    @assert hall.name in house.parameters "$(hall.name) not a known parameter"
    kk = findfirst((house.paramcomps .== hall.component) & (house.parameters .== hall.name))
    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    @assert length(hall.f) == paramspans[kk] "Length of parameter $(hall.name) unexpected: $(length(hall.f)) <> $(paramspans[kk])"
    before = sum(paramspans[1:kk-1])
    for ii in 1:paramspans[kk]
        #@assert house.b[before+ii] == 0 "Overwrite existing gradient in setobjective!"
        house.lowers[before+ii] = hall.f[ii]
    end
end

function setupper!(house::LinearProgrammingHouse, hall::LinearProgrammingHall)
    @assert hall.name in house.parameters "$(hall.name) not a known parameter"
    kk = findfirst((house.paramcomps .== hall.component) & (house.parameters .== hall.name))
    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    @assert length(hall.f) == paramspans[kk] "Length of parameter $(hall.name) unexpected: $(length(hall.f)) <> $(paramspans[kk])"
    before = sum(paramspans[1:kk-1])
    for ii in 1:paramspans[kk]
        #@assert house.b[before+ii] == 0 "Overwrite existing gradient in setobjective!"
        house.uppers[before+ii] = hall.f[ii]
    end
end

function gethouse(house::LinearProgrammingHouse, rr::Int64, cc::Int64)
    # Determine the row and column names
    varii = findlast(cumsum(varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)) .< rr) + 1
    parii = findlast(cumsum(varlengths(house.model, house.paramcomps, house.parameters)) .< cc) + 1
    rrrelative = rr - sum(varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)[1:varii-1])
    ccrelative = cc - sum(varlengths(house.model, house.paramcomps, house.parameters)[1:parii-1])
    vardims = getdims(house.model, house.constcomps[varii], get(house.constdictionary, house.constraints[varii], house.constraints[varii]))
    pardims = getdims(house.model, house.paramcomps[parii], house.parameters[parii])

    println("$(house.constcomps[varii]).$(house.constraints[varii])$(toindex(rrrelative, vardims)), $(house.paramcomps[parii]).$(house.parameters[parii])$(toindex(ccrelative, pardims)) = $(house.A[rr, cc])")
end

function constraining(house::LinearProgrammingHouse, solution::Vector{Float64}; subset=[])
    # Determine which constraint (if any) is stopping an increase or decrease of each
    if subset == []
        df = DataFrame(solution=solution)
        df[:component] = :na
        df[:parameter] = :na
    else
        df = DataFrame(solution=solution[subset])
    end
    df[:abovefail] = ""
    df[:belowfail] = ""

    # Produce names for all constraints
    varlens = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)
    names = ["" for ii in 1:sum(varlens)]
    for kk in 1:length(house.constcomps)
        ii0 = sum(varlens[1:kk-1])
        for ii in 1:varlens[kk]
            names[ii0 + ii] = "$(house.constraints[kk]).$ii"
        end
    end

    varlens = varlengths(house.model, house.paramcomps, house.parameters)
    baseconsts = house.A * solution

    println("Ignore:")
    println(join(names[find(baseconsts .> house.b)], ", "))
    ignore = baseconsts .> house.b

    if subset == []
        for kk in 1:length(house.paramcomps)
            println(kk / length(house.paramcomps))
            ii0 = sum(varlens[1:kk-1])
            for ii in 1:varlens[kk]
                df[ii0 + ii, :component] = house.paramcomps[kk]
                df[ii0 + ii, :parameter] = house.parameters[kk]

                newconst = baseconsts + house.A[:, ii0 + ii] * 1e-6 + (house.A[:, ii0 + ii] .> 0) * 1e-6
                df[ii0 + ii, :abovefail] = join(names[find((newconst .> house.b) & !ignore)], ", ")

                newconst = baseconsts - house.A[:, ii0 + ii] * 1e-6 - (house.A[:, ii0 + ii] .> 0) * 1e-6
                df[ii0 + ii, :belowfail] = join(names[find((newconst .> house.b) & !ignore)], ", ")
            end
        end
    else
        for ii in 1:length(subset)
            println(ii / length(subset))

            newconst = baseconsts + house.A[:, subset[ii]] * 1e-6 + 1e-6
            df[ii, :abovefail] = join(names[find((newconst .> house.b) & !ignore)], ", ")

            newconst = baseconsts - house.A[:, subset[ii]] * 1e-6 - 1e-6
            df[ii, :belowfail] = join(names[find((newconst .> house.b) & !ignore)], ", ")
        end
    end

    df
end

houseoptimize(house::LinearProgrammingHouse) = linprog(-house.f, house.A, '<', house.b, house.lowers, house.uppers)
houseoptimize(house::LinearProgrammingHouse, solver) = linprog(-house.f, house.A, '<', house.b, house.lowers, house.uppers, solver)
houseoptimize(house::LinearProgrammingHouse, solver, subset::Vector{Int64}) = linprog(-house.f, house.A[subset, :], '<', house.b[subset], house.lowers, house.uppers, solver)

doc"""
    findinfeasiblepair(house, solver)

Finds a range within the matrix for which the results become minimally
infeasible.  In other words, suppose that the full linear programming
matrix is $A$.  It returns $i$, $j$, such that $A[1:i, :]$ is
infeasible, but $A[1:i-1, :]$ is not, and $A[j:end, :]$ is infeasible
but $A[j+1:end, :]$ is not.

# Arguments
* `house::LinearProgrammingHouse`: An infeasible LinearProgrammingHouse.
* `solver`: A solver object which finds the infeasibility.
"""
function findinfeasiblepair(house::LinearProgrammingHouse, solver)
    top = findinfeasiblepairhelper(house, solver, "top", 2, length(house.b))
    bottom = findinfeasiblepairhelper(house, solver, "bottom", 1, length(house.b) - 1)

    [top, bottom]
end

"""
Look for the row that makes the problem feasible
"""
function findinfeasiblepairhelper(house::LinearProgrammingHouse, solver, direction, searchtop, searchbottom)
    searchpoint = round(Int, (searchtop + searchbottom) / 2 + .01) # Always round up
    println([direction, searchtop, searchbottom])
    isfeasible = checkfeasibility(house, solver, direction, searchpoint)

    if !isfeasible
        # Look for a shorter span
        if searchpoint == searchbottom
            if direction == "bottom" || checkfeasibility(house, solver, direction, searchtop)
                searchbottom
            else
                searchtop
            end
        elseif direction == "top"
            findinfeasiblepairhelper(house, solver, direction, searchtop, searchpoint)
        else
            findinfeasiblepairhelper(house, solver, direction, searchpoint, searchbottom)
        end
    else
        # Look for a longer span
        if searchpoint == searchbottom
            searchtop
        elseif (direction == "top")
            findinfeasiblepairhelper(house, solver, direction, searchpoint, searchbottom)
        else
            findinfeasiblepairhelper(house, solver, direction, searchtop, searchpoint)
        end
    end
end

function checkfeasibility(house::LinearProgrammingHouse, solver, direction, searchpoint)
    if (direction == "top")
        sol = linprog(-house.f, house.A[1:searchpoint, :], '<', house.b[1:searchpoint], house.lowers, house.uppers, solver)
    else
        sol = linprog(-house.f, house.A[searchpoint:end, :], '<', house.b[searchpoint:end], house.lowers, house.uppers, solver)
    end

    sol.status != :Infeasible
end


function summarizeparameters(house::LinearProgrammingHouse, solution::Vector{Float64})
    # Look at parameter values
    varlens = varlengths(house.model, house.paramcomps, house.parameters)
    for ii in 1:length(house.parameters)
        println(house.parameters[ii])
        index1 = sum(varlens[1:ii-1]) + 1
        index2 = sum(varlens[1:ii])

        values = solution[index1:index2]

        if (sum(values .!= 0) == 0)
            println("All zero.")
        else
            println(values[1:min(100, index2 - index1 + 1)])
            println("Sum: $(sum(values))")
        end
    end
end

"""
    getparametersolution(house, solution, parameter)

Return the array of solution values for a given parameter.
"""
function getparametersolution(house::LinearProgrammingHouse, solution::Vector{Float64}, parameter::Symbol)
    varlens = varlengths(house.model, house.paramcomps, house.parameters)

    ii = find(house.parameters .== parameter)[1]

    index1 = sum(varlens[1:ii-1]) + 1
    index2 = sum(varlens[1:ii])

    solution[index1:index2]
end

"""
    getparametersolution(house, sol, constraint)

Return the array of solution values for a given constraint.
"""
function getconstraintsolution(house, sol, constraint)
    constvalues = house.A * sol.sol

    varlens = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)

    ii = find(house.constraints .== constraint)[1]

    index1 = sum(varlens[1:ii-1]) + 1
    index2 = sum(varlens[1:ii])

    constvalues[index1:index2]
end


## Helpers

function rangeof(m::Model, name, components, names, vardictionary::Dict{Symbol, Symbol}=Dict{Symbol, Symbol}())
    varlens = varlengths(m, components, names, vardictionary)
    kk = findfirst(name .== names)
    sum(varlens[1:kk-1])+1:sum(varlens[1:kk])
end

"Translate an offset value (+1) to an index vector."
function toindex(ii::Int64, dims::Vector{Int64})
    indexes = Vector{Int64}(length(dims))
    offset = ii - 1
    for dd in 1:length(dims)
        indexes[dd] = offset % dims[dd] + 1
        offset = floor(Int64, offset / dims[dd])
    end

    return indexes
end

"Translate an index vector to an offset (+1)."
function fromindex(index::Vector{Int64}, dims::Vector{Int64})
    offset = index[end]
    for ii in length(dims)-1:-1:1
        offset = (offset - 1) * dims[ii] + index[ii]
    end

    return offset
end

"""
Return the symbols representing each of the dimensions for this variable or parameter.
"""
function getdimnames(model::Model, component::Symbol, name::Symbol)
    meta = metainfo.getallcomps()
    if name in keys(meta[(:Main, component)].parameters)
        convert(Vector{Symbol}, meta[(:Main, component)].parameters[name].dimensions)
    else
        convert(Vector{Symbol}, meta[(:Main, component)].variables[name].dimensions)
    end
end

"Return the total span occupied by each variable or parameter."
function varlengths(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, vardictionary::Dict{Symbol, Symbol}=Dict{Symbol, Symbol}())
    Int64[prod(getdims(model, components[ii], get(vardictionary, names[ii], names[ii]))) for ii in 1:length(components)]
end

## Matrix methods

"Construct a vector of values corresponding to entries in a matrix with the given dimensions, calling gen for each element."
function vectorsingle(dims::Vector{Int64}, gen)
    f = Vector{Float64}(prod(dims))
    for ii in 1:length(f)
        f[ii] = gen(toindex(ii, dims)...)
    end

    f
end

doc"""
    matrixdiagonal(dims, gen)

Creates a matrix of dimensions $\prod \text{dims}_i$.  Call the
generate function, `gen` for all indices along the diagonal.  All
combinations of indices will be called, since the "diagonal" part is
between the rows and the columns.

# Arguments
* `dims::Vector{Int64}`: The multiple dimensions collapsed into both the rows and columns.
* `gen::Function`: A function called with an argument for each dimension.

# Examples
```jldoctest
using OptiMimi

dims = [3, 2]
A = OptiMimi.matrixdiagonal(dims, (ii, jj) -> 1)
sum(A)

# output
6.0
```
"""
function matrixdiagonal(dims::Vector{Int64}, gen::Function)
    dimlen = prod(dims)
    A = spzeros(dimlen, dimlen)
    for ii in 1:dimlen
        A[ii, ii] = gen(toindex(ii, dims)...)
    end

    A
end

"""
    matrixsingle(vardims, pardims, gen)

Call the generate function for every element.  `gen` is given the
dimensions for the rows/variable followed by the dimensions for the
columns/parameters, all in order.  It should return a Float64.
"""
function matrixsingle(vardims::Vector{Int64}, pardims::Vector{Int64}, gen)
    vardimlen = prod(vardims)
    pardimlen = prod(pardims)
    A = spzeros(vardimlen, pardimlen)
    for ii in 1:vardimlen
        for jj in 1:pardimlen
            A[ii, jj] = gen(toindex(ii, vardims)..., toindex(jj, pardims)...)
        end
    end

    A
end

"""
    matrixintersect(rowdims, coldims, rowdimnames, coldimnames, gen)

Call the `gen` function with all combinations of the shared indices.
Shared dimensions must come in the same order and at the end
of the dimensions lists.

# Arguments
* `rowdims::Vector{Int64}` and `coldims::Vector{Int64}`: the number of elements in each dimension of the row and column variables.
* `rowdimnames::Vector{Symbol}` and `coldimnames::Vector{Symbol}`: The names of the dimensions; note that the last 1 or more dimensions should be the same in both of these lists.
* `gen::Function`: a function of arguments (`A`, `indices`...), where `A` is a matrix of the un-shared indices and `indices` are multiple indexes describing all shared indices.

# Examples
```jldoctest
using OptiMimi

rowdims = [3, 2]
coldims = [4, 2]
rowdimnames = [:region, :time]
coldimnames = [:person, :time]
function gen(A, tt)
    @assert tt == 1 || tt == 2
    @assert size(A) == (3, 4)
    A[2, 3] = tt
end
A = OptiMimi.matrixintersect(rowdims, coldims, rowdimnames, coldimnames, gen)
sum(A)

# output
3.0
```
"""
function matrixintersect(rowdims::Vector{Int64}, coldims::Vector{Int64}, rowdimnames::Vector{Symbol}, coldimnames::Vector{Symbol}, gen::Function)
    A = spzeros(prod(rowdims), prod(coldims))

    # Determine shared dimensions: counting from 0 = end
    numshared = 0
    while numshared < min(length(rowdims), length(coldims))
        if rowdims[end - numshared] != coldims[end - numshared]
            break
        end
        if rowdimnames[end - numshared] != coldimnames[end - numshared]
            break
        end
        numshared += 1
    end

    sharedims = rowdims[end - numshared + 1:end]
    sharedimslen = prod(sharedims)
    for ii in 1:sharedimslen
        index = toindex(ii, sharedims)
        topii = fromindex([ones(Int64, length(rowdims) - numshared); index], rowdims)
        bottomii = fromindex([rowdims[1:length(rowdims) - numshared]; index], rowdims)
        leftii = fromindex([ones(Int64, length(coldims) - numshared); index], coldims)
        rightii = fromindex([coldims[1:length(coldims) - numshared]; index], coldims)
        subA = sub(A, topii:bottomii, leftii:rightii)
        gen(subA, index...)
    end

    A
end

"""
Create an empty matrix for the given variables
"""
function matrixempty(vardims::Vector{Int64}, pardims::Vector{Int64})
    vardimlen = prod(vardims)
    pardimlen = prod(pardims)
    spzeros(vardimlen, pardimlen)
end
