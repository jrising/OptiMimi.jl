using MathProgBase
using DataFrames
using Clp

import Mimi.metainfo
import Base.*, Base.-, Base.+, Base./, Base.max

export LinearProgrammingHall, LinearProgrammingShaft, LinearProgrammingRoom, LinearProgrammingHouse
export hallsingle, hall_relabel, hallvalues, hall_duplicate
export shaftsingle, shaftvalues
export roomdiagonal, roomdiagonalintersect, roomsingle, roomchunks, roomempty, roomintersect, room_relabel, room_relabel_parameter, room_duplicate
export varsum, fromindex, fromindexes, discounted
export setobjective!, setconstraint!, setconstraintoffset!, getconstraintoffset, clearconstraint!, setlower!, setupper!, addparameter!, addconstraint!
## Debugging
export getroom, getnonzerorooms
export gethouse, constraining, houseoptimize, summarizeparameters, findinfeasiblepair, varlengths, getconstraintsolution, getparametersolution

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
function hallsingle(model::Model, component::Symbol, name::Symbol, gen::Function, dupover::Vector{Symbol}=Symbol[])
    #println("$(now()) hs $component:$name")
    if isempty(dupover)
        LinearProgrammingHall(component, name, vectorsingle(getdims(model, component, name), gen))
    else
        dupover2 = Bool[dimname in dupover for dimname in getdimnames(model, component, name)]
        LinearProgrammingHall(component, name, vectorsingle(getdims(model, component, name), gen, dupover2))
    end
end

function hallsingle(model::Model, component::Symbol, name::Symbol, constant::Float64)
    LinearProgrammingHall(component, name, vectorsingle(getdims(model, component, name), constant))
end

function hallvalues(model::Model, component::Symbol, name::Symbol, values::Array{Float64})
    #println("$(now()) hv $component:$name")
    gen(inds...) = values[inds...]
    hallsingle(model, component, name, gen)
end

"""Connect a derivative to another component: change the variable component and name to another component."""
function hall_relabel(hall::LinearProgrammingHall, from::Symbol, tocomponent::Symbol, toname::Symbol)
    @assert hall.name == from "Name mismatch in hall_relabel: $(hall.name) <> $from"

    LinearProgrammingHall(tocomponent, toname, hall.f)
end

"""Connect a derivative to another component while adding an additional dimension"""
function hall_duplicate(hall::LinearProgrammingHall, from::Symbol, tocomponent::Symbol, toname::Symbol, model::Model, newdims::Vector{Bool})
    @assert hall.name == from "Name mismatch in hall_duplicate: $(hall.name) <> $from"

    alldims = getdims(model, tocomponent, toname)

    f = zeros(prod(alldims))
    for kk in 1:prod(alldims[newdims])
        fulliis = expanddims(collect(1:length(hall.f)), alldims, newdims, toindex(kk, alldims[newdims]))
        f[fulliis] = hall.f
    end

    LinearProgrammingHall(tocomponent, toname, f)
end

function -(hall::LinearProgrammingHall)
    LinearProgrammingHall(hall.component, hall.name, -hall.f)
end

function *(hall::LinearProgrammingHall, mm::Number)
    LinearProgrammingHall(hall.component, hall.name, hall.f * mm)
end

function /(hall::LinearProgrammingHall, dd::Number)
    LinearProgrammingHall(hall.component, hall.name, hall.f / dd)
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

function shaftsingle(model::Model, component::Symbol, name::Symbol, gen::Function, dupover::Vector{Symbol}=Symbol[])
    #println("$(now()) ss $component:$name")
    if isempty(dupover)
        LinearProgrammingShaft(component, name, vectorsingle(getdims(model, component, name), gen))
    else
        dupover2 = Bool[dimname in dupover for dimname in getdimnames(model, component, name)]
        LinearProgrammingShaft(component, name, vectorsingle(getdims(model, component, name), gen, dupover2))
    end
end

function shaftsingle(model::Model, component::Symbol, name::Symbol, constant::Float64)
    LinearProgrammingShaft(component, name, vectorsingle(getdims(model, component, name), constant))
end

function shaftvalues(model::Model, component::Symbol, name::Symbol, values::Array{Float64})
    #println("$(now()) sv $component:$name")
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
function roomdiagonal(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, gen::Function, dupover::Vector{Symbol}=Symbol[])
    #println("$(now()) rd $component:$variable x $parameter")
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)

    if isempty(dupover)
        @assert dimsvar == dimspar "Variable and parameter in roomdiagonal do not have the same dimensions: $component.$variable $dimsvar <> $component.$parameter $dimspar"

        LinearProgrammingRoom(component, variable, component, parameter, matrixdiagonal(dimsvar, gen))
    else
        dupoverpar = Bool[dimname in dupover for dimname in getdimnames(model, component, parameter)]
        dupovervar = Bool[dimname in dupover for dimname in getdimnames(model, component, variable)]

        @assert dimsvar[.!dupovervar] == dimspar[.!dupoverpar] "Variable and parameter in roomdiagonal do not have the same dimensions: $component.$variable $dimsvar <> $component.$parameter $dimspar, after dup-ignore $dupover"

        matrix = matrixdiagonal(dimsvar[.!dupovervar], gen)
        matrix2 = matrixduplicate(matrix, dimsvar, dimspar, dupovervar, dupoverpar)

        LinearProgrammingRoom(component, variable, component, parameter, matrix2)
    end
end

function roomdiagonal(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, constant::Float64, dupover::Vector{Symbol}=Symbol[])
    #println("$(now()) rd $component:$variable x $parameter")
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)

    if isempty(dupover)
        @assert dimsvar == dimspar "Variable and parameter in roomdiagonal do not have the same dimensions: $component.$variable $dimsvar <> $component.$parameter $dimspar"

        LinearProgrammingRoom(component, variable, component, parameter, matrixdiagonal(dimsvar, constant))
    else
        dupoverpar = Bool[dimname in dupover for dimname in getdimnames(model, component, parameter)]
        dupovervar = Bool[dimname in dupover for dimname in getdimnames(model, component, variable)]

        @assert dimsvar[.!dupovervar] == dimspar[.!dupoverpar] "Variable and parameter in roomdiagonal do not have the same dimensions: $component.$variable $dimsvar <> $component.$parameter $dimspar, after dup-ignore $dupover"

        matrix = matrixdiagonal(dimsvar[.!dupovervar], constant)
        matrix2 = matrixduplicate(matrix, dimsvar, dimspar, dupovervar, dupoverpar)
        LinearProgrammingRoom(component, variable, component, parameter, matrix2)
    end
end

"""
gen returns a full matrix for all shared variables, in variable dimension order.
"""
function roomdiagonalintersect(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, gen::Function)
    #println("$(now()) rd $component:$variable x $parameter")
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)
    dimnamesvar = getdimnames(model, component, variable)
    dimnamespar = getdimnames(model, component, parameter)

    A = matrixdiagonalintersect(dimsvar, dimspar, dimnamesvar, dimnamespar, gen)
    LinearProgrammingRoom(component, variable, component, parameter, A)
end

"""
Fill in every element
"""
function roomsingle(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, gen::Function, vardupover::Vector{Symbol}=Symbol[], pardupover::Vector{Symbol}=Symbol[])
    #println("$(now()) rs $component:$variable x $parameter")
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)

    if isempty(vardupover) && isempty(pardupover)
        LinearProgrammingRoom(component, variable, component, parameter, matrixsingle(dimsvar, dimspar, gen))
    else
        vardupover2 = Bool[dimname in vardupover for dimname in getdimnames(model, component, variable)]
        pardupover2 = Bool[dimname in pardupover for dimname in getdimnames(model, component, parameter)]
        LinearProgrammingRoom(component, variable, component, parameter, matrixsingle(dimsvar, dimspar, gen, vardupover2, pardupover2))
    end
end

"""
Call with matrices for each combination of chunk variables
    Chunked variables must be consecutive at the end.
"""
function roomchunks(model::Model, component::Symbol, variable::Symbol, parameter::Symbol, gen::Function, varchunk::Vector{Symbol}, parchunk::Vector{Symbol})
    #println("$(now()) rc $component:$variable x $parameter")
    dimsvar = getdims(model, component, variable)
    dimspar = getdims(model, component, parameter)

    dimvarnames = getdimnames(model, component, variable)
    dimparnames = getdimnames(model, component, parameter)

    varchunk2 = 0
    while varchunk2 < length(dimsvar)
        if !(dimvarnames[end - varchunk2] in varchunk)
            break
        end
        varchunk2 += 1
    end

    parchunk2 = 0
    while parchunk2 < length(dimspar)
        if !(dimparnames[end - parchunk2] in parchunk)
            break
        end
        parchunk2 += 1
    end

    LinearProgrammingRoom(component, variable, component, parameter, matrixchunks(dimsvar, dimspar, gen, varchunk2, parchunk2))
end

"""
    roomintersect(model, component, variable1, variable2, gen)

Generate a room at the intersection of two variables.  Variable1 is
the constraint, and variable2 has the optimization parameter.

Call `gen` for each shared index, passing an array to be filled for unshared indexes.
This version assumes both variables come from the same component.

See `matrixintersect` for the matrix generation logic.
"""
function roomintersect(model::Model, component::Symbol, variable1::Symbol, variable2::Symbol, gen::Function, dupover1::Vector{Symbol}=Symbol[], dupover2::Vector{Symbol}=Symbol[])
    #println("$(now()) ri $component:$variable1 x $variable2")
    dims1 = getdims(model, component, variable1)
    dims2 = getdims(model, component, variable2)
    dimnames1 = getdimnames(model, component, variable1)
    dimnames2 = getdimnames(model, component, variable2)

    if isempty(dupover1) && isempty(dupover2)
        LinearProgrammingRoom(component, variable1, component, variable2, matrixintersect(dims1, dims2, dimnames1, dimnames2, gen))
    else
        dupover12 = Bool[dimname in dupover1 for dimname in getdimnames(model, component, variable1)]
        dupover22 = Bool[dimname in dupover2 for dimname in getdimnames(model, component, variable2)]
        LinearProgrammingRoom(component, variable1, component, variable2, matrixintersect(dims1, dims2, dimnames1, dimnames2, gen, dupover12, dupover22))
    end
end

"""
    roomintersect(model, component, variable1, component2, variable2, gen)

Generate a room at the intersection of two variables.

Call `gen` for each shared index, passing an array to be filled for unshared indexes.
This version allows the variables to come from different component.

See `matrixintersect` for the matrix generation logic.
"""
function roomintersect(model::Model, component1::Symbol, variable1::Symbol, component2::Symbol, variable2::Symbol, gen::Function, dupover1::Vector{Symbol}=Symbol[], dupover2::Vector{Symbol}=Symbol[])
    #println("$(now()) ri $component1:$variable1 x $component2:$variable2")
    dims1 = getdims(model, component1, variable1)
    dims2 = getdims(model, component2, variable2)
    dimnames1 = getdimnames(model, component1, variable1)
    dimnames2 = getdimnames(model, component2, variable2)

    if isempty(dupover1) && isempty(dupover2)
        LinearProgrammingRoom(component1, variable1, component2, variable2, matrixintersect(dims1, dims2, dimnames1, dimnames2, gen))
    else
        dupover12 = Bool[dimname in dupover1 for dimname in getdimnames(model, component1, variable1)]
        dupover22 = Bool[dimname in dupover2 for dimname in getdimnames(model, component2, variable2)]
        LinearProgrammingRoom(component1, variable1, component2, variable2, matrixintersect(dims1, dims2, dimnames1, dimnames2, gen, dupover12, dupover22))
    end
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

"""Connect a derivative to another component while adding an additional dimension"""
function room_duplicate(room::LinearProgrammingRoom, tovariable::Symbol, toparameter::Symbol, model::Model)
    parolddimnames = getdimnames(model, room.paramcomponent, room.parameter)
    parnewdimnames = getdimnames(model, room.paramcomponent, toparameter)
    @assert issubset(parolddimnames, parnewdimnames)
    varolddimnames = getdimnames(model, room.varcomponent, room.variable)
    varnewdimnames = getdimnames(model, room.varcomponent, tovariable)
    @assert issubset(varolddimnames, varnewdimnames)

    rowdupover = Bool[]
    jj = 1
    for ii in 1:length(varnewdimnames)
        if jj <= length(varolddimnames) && varnewdimnames[ii] == varolddimnames[jj]
            jj += 1
            push!(rowdupover, false)
        else
            push!(rowdupover, true)
        end
    end
    @assert jj == length(varolddimnames) + 1 "Got $jj after matching $varolddimnames in $varnewdimnames for $rowdupover"

    coldupover = Bool[]
    jj = 1
    for ii in 1:length(parnewdimnames)
        if jj <= length(parolddimnames) && parnewdimnames[ii] == parolddimnames[jj]
            jj += 1
            push!(coldupover, false)
        else
            push!(coldupover, true)
        end
    end
    @assert jj == length(parolddimnames) + 1 "Got $jj after matching $parolddimnames in $parnewdimnames for $coldupover"

    A = matrixduplicate(room.A, getdims(model, room.varcomponent, tovariable), getdims(model, room.paramcomponent, toparameter), rowdupover, coldupover)
    LinearProgrammingRoom(room.varcomponent, tovariable, room.paramcomponent, toparameter, A)
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

function varsum(room::LinearProgrammingRoom, axis::Int64, model::Model, newvar::Symbol)
    vardims = getdims(model, room.varcomponent, room.variable)
    iis, jjs, vvs = findnz(room.A)

    indexes = toindexes(iis, vardims)
    newiis = fromindexes(indexes[:, 1:length(vardims) .!= axis], vardims[1:length(vardims) .!= axis])

    useiis = []
    usejjs = []
    usevvs = []
    for newii in unique(newiis)
        withinii = newiis .== newii
        for newjj in unique(jjs[withinii])
            withinij = jjs[withinii] .== newjj
            push!(useiis, newii)
            push!(usejjs, newjj)
            push!(usevvs, sum(vvs[withinii][withinij]))
        end
    end

    A = sparse(useiis, usejjs, usevvs, prod(vardims[1:length(vardims) .!= axis]), size(room.A, 2))
    LinearProgrammingRoom(room.varcomponent, newvar, room.paramcomponent, room.parameter, A)
end

"""
Discount economic costs in room.  Only rows are discounted (i.e., doesn't matter when action made, just when cost incurred).
"""
function discounted(model::Model, room::LinearProgrammingRoom, rate::Float64)
    dims = getdims(model, room.varcomponent, room.variable)

    dupover2 = Bool[dimname != :time for dimname in getdimnames(model, room.varcomponent, room.variable)]
    outers, inners = interpretdupover(dupover2)
    dupouter = prod(dims[outers])
    dupinner = prod(dims[inners])

    timevalues = collect(1:getindexcount(model, :time))
    alltime = repeat(timevalues, inner=[dupinner], outer=[dupouter])

    discount = exp.(-rate * (alltime - 1))
    # discounted = room.A .* discount
    iis, jjs, vvs = findnz(room.A)
    vvs2 = vvs .* discount[iis]
    discounted = sparse(iis, jjs, vvs2, size(room.A, 1), size(room.A, 2))

    LinearProgrammingRoom(room.varcomponent, room.variable, room.paramcomponent, room.parameter, discounted)
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
    kk = findfirst((house.paramcomps .== hall.component) .& (house.parameters .== hall.name))
    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    @assert length(hall.f) == paramspans[kk] "Length of parameter $(hall.name) unexpected: $(length(hall.f)) <> $(paramspans[kk])"
    before = sum(paramspans[1:kk-1])
    for ii in 1:paramspans[kk]
        #@assert house.f[before+ii] == 0 "Overwrite existing gradient in setobjective!"
        house.f[before+ii] = hall.f[ii]
    end
end

function setconstraint!(house::LinearProgrammingHouse, room::LinearProgrammingRoom)
    kk = findfirst((house.paramcomps .== room.paramcomponent) .& (house.parameters .== room.parameter))
    @assert kk > 0 "$(room.paramcomponent).$(room.parameter) not a known parameter"
    ll = findfirst((house.constcomps .== room.varcomponent) .& (house.constraints .== room.variable))
    @assert ll > 0 "$(room.varcomponent).$(room.variable) not a known variable"

    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)
    @assert size(room.A, 1) == constspans[ll] "Length of variable $(room.variable) unexpected: $(size(room.A, 1)) given <> $(constspans[ll])"
    @assert size(room.A, 2) == paramspans[kk] "Length of parameter $(room.parameter) unexpected: $(size(room.A, 2)) given <> $(paramspans[kk])"

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
    rowkk = findfirst((house.constcomps .== varcomponent) .& (house.constraints .== variable))
    rowbefore = sum(constspans[1:rowkk-1])

    paramspans = varlengths(house.model, house.paramcomps, house.parameters)
    colkk = findfirst((house.paramcomps .== paramcomponent) .& (house.parameters .== parameter))
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
    kk = findfirst((house.constcomps .== component) .& (house.constraints .== variable))
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
    kk = findfirst((house.constcomps .== component) .& (house.constraints .== variable))
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

    ll = findfirst((house.constcomps .== component) .& (house.constraints .== variable))
    constspans = varlengths(house.model, house.constcomps, house.constraints, house.constdictionary)

    constbefore = sum(constspans[1:ll-1])

    house.A[constbefore+1:constbefore+constspans[ll], :] = 0
    house.b[constbefore+1:constbefore+constspans[ll]] = 0
end

function setlower!(house::LinearProgrammingHouse, hall::LinearProgrammingHall)
    @assert hall.name in house.parameters "$(hall.name) not a known parameter"
    kk = findfirst((house.paramcomps .== hall.component) .& (house.parameters .== hall.name))
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
    kk = findfirst((house.paramcomps .== hall.component) .& (house.parameters .== hall.name))
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
                df[ii0 + ii, :abovefail] = join(names[find((newconst .> house.b) .& .!ignore)], ", ")

                newconst = baseconsts - house.A[:, ii0 + ii] * 1e-6 - (house.A[:, ii0 + ii] .> 0) * 1e-6
                df[ii0 + ii, :belowfail] = join(names[find((newconst .> house.b) .& .!ignore)], ", ")
            end
        end
    else
        for ii in 1:length(subset)
            println(ii / length(subset))

            newconst = baseconsts + house.A[:, subset[ii]] * 1e-6 + 1e-6
            df[ii, :abovefail] = join(names[find((newconst .> house.b) .& !ignore)], ", ")

            newconst = baseconsts - house.A[:, subset[ii]] * 1e-6 - 1e-6
            df[ii, :belowfail] = join(names[find((newconst .> house.b) .& !ignore)], ", ")
        end
    end

    df
end

houseoptimize(house::LinearProgrammingHouse) = linprog(-house.f, house.A, '<', house.b, house.lowers, house.uppers, ClpSolver())
houseoptimize(house::LinearProgrammingHouse, solver) = linprog(-house.f, house.A, '<', house.b, house.lowers, house.uppers, solver)
houseoptimize(house::LinearProgrammingHouse, solver, subset::Vector{Int64}) = linprog(-house.f, house.A[subset, :], '<', house.b[subset], house.lowers, house.uppers, solver)

doc"""
    findinfeasiblepair(house, solver)

Finds a range within the matrix for which the results become minimally
infeasible.  In other words, suppose that the full linear programming
matrix is $A$.  It returns $i$, $j$, such that $A[1:i, :]$ is
infeasible, but $A[1:i-1, :]$ is not, and $A[j:end, :]$ is infeasible
but $A[j-1:end, :]$ is not.

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

"""
Translate an offset value (+1) to an index vector.
This faster if dims is a vector.  Use ind2sub if dims is tuple.
"""
function toindex(ii::Int64, dims::Vector{Int64})
    indexes = Vector{Int64}(length(dims))
    offset = ii - 1
    for dd in 1:length(dims)
        indexes[dd] = offset % dims[dd] + 1
        offset = div(offset, dims[dd])
    end

    return indexes
end

"""
As `toindex`, but for many iis.  Ruturns N x D
"""
function toindexes(iis::Vector{Int64}, dims::Vector{Int64})
    indexes = Matrix{Int64}(length(iis), length(dims))
    offsets = iis - 1
    for dd in 1:length(dims)
        indexes[:, dd] = offsets .% dims[dd] + 1
        offsets = div.(offsets, dims[dd])
    end

    return indexes
end

"""
Translate an index vector to an offset (+1).
This seems to be faster than sub2ind irrespective of the tupleness of dims.
"""
function fromindex(index::Vector{Int64}, dims::Vector{Int64})
    offset = index[end]
    for ii in length(dims)-1:-1:1
        offset = (offset - 1) * dims[ii] + index[ii]
    end

    return offset
end

"""
As `fromindex` but for many indexes (N x D)
"""
function fromindexes(indexes::Matrix{Int64}, dims::Vector{Int64})
    offsets = indexes[:, end]
    for ii in length(dims)-1:-1:1
        offsets = (offsets - 1) * dims[ii] + indexes[:, ii]
    end

    return offsets
end

"""
Insert an extra dimension somewhere, adjusting offset (+1) values accordingly
"""
function insertdim(iis::Vector{Int64}, dims::Vector{Int64}, insertat::Int64, dimvalue::Int64, dimsize::Int64)
    @assert insertat >= 1 && insertat <= length(dims) + 1
    @assert dimvalue >= 1 && dimvalue <= dimsize
    if insertat == 1
        (iis - 1) * dimsize + dimvalue
    elseif insertat == length(dims) + 1
        (dimvalue - 1) * prod(dims) + iis
    else
        proddimsleft = prod(dims[1:insertat-1])
        (div.(iis - 1, proddimsleft) * dimsize + dimvalue - 1) * proddimsleft + (iis - 1) .% proddimsleft + 1
    end
end

"""
Insert as many dims as needed to expand to the alldims
"""
function expanddims(iis::Vector{Int64}, alldims::Vector{Int64}, newdims::AbstractVector{Bool}, newdimvalues::Vector{Int64})
    insertat = findfirst(newdims)
    iis2 = insertdim(iis, alldims[.!newdims], insertat, newdimvalues[1], alldims[insertat])
    if length(newdimvalues) == 1
        return iis2
    else
        newdims[insertat] = false
        return expanddims(iis2, alldims, newdims, newdimvalues[2:end])
    end
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

"Splits dupover into inner and outer dimensions"
function interpretdupover(dupover::Vector{Bool})
    ii = 1

    # Initial trues
    inners = Int64[]
    if dupover[1]
        while ii < length(dupover)
            push!(inners, ii)
            ii += 1
            if !dupover[ii]
                break
            end
        end
    end
    if dupover[ii] # fell off the end
        error("Cannot duplicate all dimensions.")
    end

    # Middle falses
    while ii < length(dupover)
        ii += 1
        if dupover[ii]
            break
        end
    end
    if !dupover[ii]
        return Int64[], inners
    end

    # Final trues
    outers = Int64[]
    while ii <= length(dupover)
        if !dupover[ii]
            error("Cannot support middle dimension duplication.")
        end

        push!(outers, ii)
        ii += 1
    end

    return outers, inners
end

"Construct a vector of values corresponding to entries in a matrix with the given dimensions, calling gen for each element."
function vectorsingle(dims::Vector{Int64}, gen::Function)
    f = Vector{Float64}(prod(dims))
    for ii in 1:length(f)
        f[ii] = gen(toindex(ii, dims)...)
    end

    f
end

function vectorsingle(dims::Vector{Int64}, constant::Float64)
    return constant * ones(prod(dims))
end

function vectorsingle(dims::Vector{Int64}, gen::Function, dupover::Vector{Bool})
    outers, inners = interpretdupover(dupover)
    dupouter = prod(dims[outers])
    dupinner = prod(dims[inners])

    f = Vector{Float64}(prod(dims[.!dupover]))
    for ii in 1:length(f)
        f[ii] = gen(toindex(ii, dims[.!dupover])...)
    end

    repeat(f, inner=[dupinner], outer=[dupouter])
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
    index = ones(Int64, length(dims))
    for ii in 1:dimlen
        A[ii, ii] = gen(index...)

        index[1] += 1
        for kk in 1:(length(dims) - 1)
            if index[kk] <= dims[kk]
                break
            end
            index[kk] = 1
            index[kk+1] += 1
        end
        #A[ii, ii] = gen(toindex(ii, dims)...) # TOO SLOW
    end

    A
end

function matrixdiagonal(dims::Vector{Int64}, constant::Float64)
    dimlen = prod(dims)
    speye(dimlen, dimlen) * constant
end

function matrixdiagonal(dims::Vector{Int64}, gen::Function, dupover::Vector{Bool})
    diag = vectorsingle(dims, gen, dupover)

    dimlen = prod(dims)
    A = spzeros(dimlen, dimlen)
    for ii in 1:dimlen
        if diag[ii] != 0
            A[ii, ii] = diag[ii]
        end
    end

    A
end

"""
Identifies all the common dimensions, calling the generator function
with an intersection matrix once for each value of the shared
dimensions (filling out their diagonal).

gen is called with each combination of unshared variables, and
returns a full matrix for all shared variables, in variable dimension order.
"""
function matrixdiagonalintersect(rowdims::Vector{Int64}, coldims::Vector{Int64}, rowdimnames::Vector{Symbol}, coldimnames::Vector{Symbol}, gen::Function)
    common = intersect(rowdimnames, coldimnames)

    commonrow = Bool[dimname in common for dimname in rowdimnames]
    commoncol = Bool[dimname in common for dimname in coldimnames]

    remainrowdims = rowdims[.!commonrow]
    remaincoldims = coldims[.!commoncol]

    fullrowsub = repmat([0], length(rowdims))
    fullcolsub = repmat([0], length(coldims))

    alliis = []
    alljjs = []
    allvvs = []
    for ii in 1:prod(remainrowdims)
        for jj in 1:prod(remaincoldims)
            remainrowsub = toindex(ii, remainrowdims)
            remaincolsub = toindex(jj, remaincoldims)

            Apart = gen(remainrowsub..., remaincolsub...)
            kks = collect(1:prod(rowdims[commonrow]))

            if (sum(.!commonrow) == 0)
                append!(alliis, kks)
            else
                append!(alliis, expanddims(kks, rowdims, .!commonrow, remainrowsub))
            end
            if (sum(.!commoncol) == 0)
                append!(alljjs, kks)
            else
                append!(alljjs, expanddims(kks, coldims, .!commoncol, remaincolsub))
            end
            append!(allvvs, vec(Apart))
        end
    end

    sparse(alliis, alljjs, allvvs, prod(rowdims), prod(coldims))
end


"""
    matrixsingle(vardims, pardims, gen)

Call the generate function for every element.  `gen` is given the
dimensions for the rows/variable followed by the dimensions for the
columns/parameters, all in order.  It should return a Float64.
"""
function matrixsingle(vardims::Vector{Int64}, pardims::Vector{Int64}, gen::Function)
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

function matrixsingle(vardims::Vector{Int64}, pardims::Vector{Int64}, gen::Function, vardupover::Vector{Bool}, pardupover::Vector{Bool})
    outers, inners = interpretdupover(vardupover)
    vardupouter = prod(vardims[outers])
    vardupinner = prod(vardims[inners])

    outers, inners = interpretdupover(pardupover)
    pardupouter = prod(pardims[outers])
    pardupinner = prod(pardims[inners])

    vardimlen = prod(vardims[.!vardupover])
    pardimlen = prod(pardims[.!pardupover])
    A = zeros(vardimlen, pardimlen) # Not sparse
    for ii in 1:vardimlen
        for jj in 1:pardimlen
            A[ii, jj] = gen(toindex(ii, vardims[.!vardupover])..., toindex(jj, pardims[.!pardupover])...)
        end
    end

    sparse(repeat(A, inner=[vardupinner, pardupinner], outer=[vardupouter, pardupouter]))
end

"""
matrixchunks(rowdims, coldims, gen, rowchunks, colchunks)
Call the `gen` function with all combinations of chunked indices.
    rowchunks and colchunks is number of row and col vars at end to chunk.
"""
function matrixchunks(rowdims::Vector{Int64}, coldims::Vector{Int64}, gen::Function, rowchunk::Int64, colchunk::Int64)
    alliis = Int64[]
    alljjs = Int64[]
    allvvs = Float64[]
    #A = spzeros(prod(rowdims), prod(coldims))

    chunkrowdims = rowdims[end-rowchunk+1:end]
    chunkrowdimlen = prod(chunkrowdims)
    chunkcoldims = coldims[end-colchunk+1:end]
    chunkcoldimlen = prod(chunkcoldims)
    for ii in 1:chunkrowdimlen
        rowindex = toindex(ii, chunkrowdims)
        for jj in 1:chunkcoldimlen
            colindex = toindex(jj, chunkcoldims)

            topii = fromindex([ones(Int64, length(rowdims) - rowchunk); rowindex], rowdims)
            bottomii = fromindex([rowdims[1:length(rowdims) - rowchunk]; rowindex], rowdims)
            leftii = fromindex([ones(Int64, length(coldims) - colchunk); colindex], coldims)
            rightii = fromindex([coldims[1:length(coldims) - colchunk]; colindex], coldims)

            subA = gen(rowindex..., colindex...)
            @assert size(subA)[1] == bottomii - topii + 1 "Got matrix of size $(size(subA)) for rows $topii - $bottomii"
            @assert size(subA)[2] == rightii - leftii + 1 "Got matrix of size $(size(subA)) for columns $leftii - $rightii"

            iis, jjs, vvs = findnz(subA)
            append!(alliis, iis + topii - 1)
            append!(alljjs, jjs + leftii - 1)
            append!(allvvs, vvs)
            #A[topii:bottomii, leftii:rightii] = gen(rowindex..., colindex...)
        end
    end

    sparse(alliis, alljjs, allvvs, prod(rowdims), prod(coldims))
    #A
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
        subA = view(A, topii:bottomii, leftii:rightii)
        gen(subA, index...)
    end

    A
end

function matrixintersect(rowdims::Vector{Int64}, coldims::Vector{Int64}, rowdimnames::Vector{Symbol}, coldimnames::Vector{Symbol}, gen::Function, rowdupover::Vector{Bool}, coldupover::Vector{Bool})
    # Figure out how many of dupovers are shared
    rowdupovershared = zeros(Bool, length(rowdupover))
    coldupovershared = zeros(Bool, length(coldupover))
    for ii in 0:(min(length(rowdims), length(coldims))-1)
        if rowdupover[end-ii] && coldupover[end-ii] && rowdimnames[end-ii] == coldimnames[end-ii]
            rowdupovershared[end-ii] = true
            coldupovershared[end-ii] = true
        else
            break
        end
    end

    # Generate, without any dups
    A = matrixintersect(rowdims[.!rowdupover], coldims[.!coldupover], rowdimnames[.!rowdupover], coldimnames[.!coldupover], gen)

    A2 = matrixduplicate(A, rowdims[.!rowdupovershared], coldims[.!coldupovershared], rowdupover[.!rowdupovershared], coldupover[.!coldupovershared])

    alliis, alljjs, allvvs = findnz(A2)

    # Create the shared-dupped portion

    shareddupouter = prod(rowdims[rowdupovershared])
    rownotshared = prod(rowdims[.!rowdupovershared])
    colnotshared = prod(coldims[.!coldupovershared])

    allvvs2 = repeat(allvvs, outer=[shareddupouter])
    alliis2 = zeros(Int64, length(alliis) * shareddupouter)
    alljjs2 = zeros(Int64, length(alljjs) * shareddupouter)

    for kk in 1:shareddupouter
        alliis2[(kk - 1) * length(alliis) + (1:length(alliis))] = (kk - 1) * rownotshared + alliis
        alljjs2[(kk - 1) * length(alljjs) + (1:length(alljjs))] = (kk - 1) * colnotshared + alljjs
    end

    sparse(alliis2, alljjs2, allvvs2, prod(rowdims), prod(coldims))
end

"""
Create an empty matrix for the given variables
"""
function matrixempty(vardims::Vector{Int64}, pardims::Vector{Int64})
    vardimlen = prod(vardims)
    pardimlen = prod(pardims)
    spzeros(vardimlen, pardimlen)
end

"""
General duplication of a known matrix over additional dimensions.
Falls back on matrixduplicate_general if not all inner-most and outer-most dupovers
"""
function matrixduplicate(A::SparseMatrixCSC{Float64, Int64}, rowdims::Vector{Int64}, coldims::Vector{Int64}, rowdupover::Vector{Bool}, coldupover::Vector{Bool})
    if sum(abs.(diff([true; rowdupover]))) <= 2 && sum(abs.(diff([true; coldupover]))) <= 2
        matrixduplicate_extremes(A, rowdims, coldims, rowdupover, coldupover)
    else
        matrixduplicate_general(A, rowdims, coldims, rowdupover, coldupover)
    end
end

"""
Duplication of a known matrix over inner-most and outer-most dimensions
"""
function matrixduplicate_extremes(A::SparseMatrixCSC{Float64, Int64}, rowdims::Vector{Int64}, coldims::Vector{Int64}, rowdupover::Vector{Bool}, coldupover::Vector{Bool})
    outers, inners = interpretdupover(rowdupover)
    rowdupouter = prod(rowdims[outers])
    rowdupinner = prod(rowdims[inners])

    outers, inners = interpretdupover(coldupover)
    coldupouter = prod(coldims[outers])
    coldupinner = prod(coldims[inners])

    # Create fully-dupped portion
    iis, jjs, vvs = findnz(A)

    kkdupnum = rowdupouter*rowdupinner*coldupouter*coldupinner

    allvvs = repeat(vvs, inner=[kkdupnum])
    # Need to fill these next two in
    alliis = zeros(Int64, length(iis) * kkdupnum)
    alljjs = zeros(Int64, length(jjs) * kkdupnum)

    nrows = size(A)[1]
    ncols = size(A)[2]

    for kk in 1:length(iis)
        ## order is vec([(slow, fast) for fast in 1:N, slow in 1:M])
        rowduped = vec([(rowoo - 1) * nrows * rowdupinner + (iis[kk] - 1) * rowdupinner + rowii for rowii in 1:rowdupinner, rowoo in 1:rowdupouter])
        colduped = vec([(coloo - 1) * ncols * coldupinner + (jjs[kk] - 1) * coldupinner + colii for colii in 1:coldupinner, coloo in 1:coldupouter])
        alliis[(kk - 1) * kkdupnum + (1:kkdupnum)] = repeat(rowduped, outer=[coldupouter*coldupinner])
        alljjs[(kk - 1) * kkdupnum + (1:kkdupnum)] = repeat(colduped, inner=[rowdupouter*rowdupinner])
    end

    sparse(alliis, alljjs, allvvs, prod(rowdims), prod(coldims))
end

"""
General duplication of a known matrix over additional dimensions.
This is not constrained to have dupover only on inner-most and outer-most dimensions.
"""
function matrixduplicate_general(origA::SparseMatrixCSC{Float64, Int64}, rowdims::Vector{Int64}, coldims::Vector{Int64}, rowdupover::Vector{Bool}, coldupover::Vector{Bool})
    iis, jjs, vvs = findnz(origA)

    alliis = Int64[]
    alljjs = Int64[]
    allvvs = Float64[]

    if !any(rowdupover)
        newdimsizes = coldims[coldupover]
        for kk in 1:prod(coldims[coldupover])
            fulljjs = expanddims(jjs, coldims, coldupover, toindex(kk, newdimsizes))
            append!(alliis, iis)
            append!(alljjs, fulljjs)
            append!(allvvs, vvs)
        end
    elseif !any(coldupover)
        newdimsizes = rowdims[rowdupover]
        for kk in 1:prod(rowdims[rowdupover])
            fulliis = expanddims(iis, rowdims, rowdupover, toindex(kk, newdimsizes))
            append!(alliis, fulliis)
            append!(alljjs, jjs)
            append!(allvvs, vvs)
        end
    else
        newrowsizes = rowdims[rowdupover]
        newcolsizes = coldims[coldupover]
        for rowkk in 1:prod(rowdims[rowdupover])
            for colkk in 1:prod(coldims[coldupover])
                fulliis = expanddims(iis, rowdims, rowdupover, toindex(rowkk, newrowsizes))
                fulljjs = expanddims(jjs, coldims, coldupover, toindex(colkk, newcolsizes))
                append!(alliis, fulliis)
                append!(alljjs, fulljjs)
                append!(allvvs, vvs)
            end
        end
    end

    sparse(alliis, alljjs, allvvs, prod(rowdims), prod(coldims))
end
