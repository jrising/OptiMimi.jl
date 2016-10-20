
<a id='Linear-Programming-Guide-1'></a>

# Linear Programming Guide


<a id='Types-1'></a>

## Types

<a id='OptiMimi.LinearProgrammingHall' href='#OptiMimi.LinearProgrammingHall'>#</a>
**`OptiMimi.LinearProgrammingHall`** &mdash; *Type*.



```
LinearProgrammingHall
```

A vector of values, used either to describes the objective function or the constraint offsets.

**Fields**

  * `component::Symbol`: The component defining a parameter or variable.
  * `name::Symbol`: Either a parameter or variable name.
  * `f::Vector{Float64}`: A vector of values, for every entry in the variable.

<a id='OptiMimi.LinearProgrammingHouse' href='#OptiMimi.LinearProgrammingHouse'>#</a>
**`OptiMimi.LinearProgrammingHouse`** &mdash; *Type*.



```
LinearProgrammingHouse
```

The full description of a linear programming problem, including all its variables, parameters, the constraints, and the objective.

The linear programming that is solved is always: $\max f' x$ $A x \le b$ $x_{lower} \le x \le x_{upper}$

**Fields**

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

<a id='OptiMimi.LinearProgrammingRoom' href='#OptiMimi.LinearProgrammingRoom'>#</a>
**`OptiMimi.LinearProgrammingRoom`** &mdash; *Type*.



```
LinearProgrammingRoom
```

A matrix of values, used to describe a part of the linear programming problem constraint matrix (a gradient matrix).  The rows of the matrix correspond to variable indices, and the columns correspond to parameter indices.  The values are stored as a sparse matrix, since most parameter indices are assumed to not interact with most variable indices.

**Fields**

  * `varcomponent::Symbol`: The component defining the variable.
  * `variable::Symbol`: The variable name.
  * `paramcomponent::Symbol`: The component defining the parameter.
  * `parameter::Symbol`: The parameter name.
  * `A::SparseMatrixCSC{Float64, Int64}`: A sparse matrix of values, for every combination of the parameter and variable.

<a id='OptiMimi.LinearProgrammingShaft' href='#OptiMimi.LinearProgrammingShaft'>#</a>
**`OptiMimi.LinearProgrammingShaft`** &mdash; *Type*.



```
LinearProgrammingShaft
```

A vector of values, used either to describes the objective function or the constraint offsets.  This acts as the transpose of a Hall, although this distinction is only important for multiplying a Room by Shaft, which then returns a Hall.

**Fields**

  * `component::Symbol`: The component defining a parameter or variable.
  * `name::Symbol`: Either a parameter or variable name.
  * `x::Vector{Float64}`: A vector of values, for every entry in the variable.


<a id='Functions-1'></a>

## Functions

<a id='OptiMimi.findinfeasiblepair-Tuple{OptiMimi.LinearProgrammingHouse,Any}' href='#OptiMimi.findinfeasiblepair-Tuple{OptiMimi.LinearProgrammingHouse,Any}'>#</a>
**`OptiMimi.findinfeasiblepair`** &mdash; *Method*.



```
findinfeasiblepair(house, solver)
```

Finds a range within the matrix for which the results become minimally infeasible.  In other words, suppose that the full linear programming matrix is $A$.  It returns $i$, $j$, such that $A[1:i, :]$ is infeasible, but $A[1:i-1, :]$ is not, and $A[j:end, :]$ is infeasible but $A[j+1:end, :]$ is not.

**Arguments**

  * `house::LinearProgrammingHouse`: An infeasible LinearProgrammingHouse.
  * `solver`: A solver object which finds the infeasibility.

<a id='OptiMimi.fromindex-Tuple{Array{Int64,1},Array{Int64,1}}' href='#OptiMimi.fromindex-Tuple{Array{Int64,1},Array{Int64,1}}'>#</a>
**`OptiMimi.fromindex`** &mdash; *Method*.



Translate an index vector to an offset (+1).

<a id='OptiMimi.getconstraintoffset-Tuple{OptiMimi.LinearProgrammingHouse,Symbol,Symbol}' href='#OptiMimi.getconstraintoffset-Tuple{OptiMimi.LinearProgrammingHouse,Symbol,Symbol}'>#</a>
**`OptiMimi.getconstraintoffset`** &mdash; *Method*.



```
getconstraintoffset(house, component, variable, reshp)
```

Return the values for a constraint, optionally reshaped to the original dimensions.

**Arguments**

  * `house::LinearProgrammingHouse`: The house from which to get the values.
  * `component::Symbol`: The component for the constraint variable.
  * `variable::Symbol`: The variable for the constraint.
  * `reshp::Bool`: Should it be reshaped to the original variable dimensions? (default: false)

<a id='OptiMimi.getconstraintsolution-Tuple{Any,Any,Any}' href='#OptiMimi.getconstraintsolution-Tuple{Any,Any,Any}'>#</a>
**`OptiMimi.getconstraintsolution`** &mdash; *Method*.



```
getparametersolution(house, sol, constraint)
```

Return the array of solution values for a given constraint.

<a id='OptiMimi.getparametersolution-Tuple{OptiMimi.LinearProgrammingHouse,Array{Float64,1},Symbol}' href='#OptiMimi.getparametersolution-Tuple{OptiMimi.LinearProgrammingHouse,Array{Float64,1},Symbol}'>#</a>
**`OptiMimi.getparametersolution`** &mdash; *Method*.



```
getparametersolution(house, solution, parameter)
```

Return the array of solution values for a given parameter.

<a id='OptiMimi.hall_relabel-Tuple{OptiMimi.LinearProgrammingHall,Symbol,Symbol,Symbol}' href='#OptiMimi.hall_relabel-Tuple{OptiMimi.LinearProgrammingHall,Symbol,Symbol,Symbol}'>#</a>
**`OptiMimi.hall_relabel`** &mdash; *Method*.



Connect a derivative to another component: change the variable component and name to another component.

<a id='OptiMimi.room_relabel-Tuple{OptiMimi.LinearProgrammingRoom,Symbol,Symbol,Symbol}' href='#OptiMimi.room_relabel-Tuple{OptiMimi.LinearProgrammingRoom,Symbol,Symbol,Symbol}'>#</a>
**`OptiMimi.room_relabel`** &mdash; *Method*.



Connect a gradient to another component: change the variable component and name to another component.

<a id='OptiMimi.room_relabel_parameter-Tuple{OptiMimi.LinearProgrammingRoom,Symbol,Symbol,Symbol}' href='#OptiMimi.room_relabel_parameter-Tuple{OptiMimi.LinearProgrammingRoom,Symbol,Symbol,Symbol}'>#</a>
**`OptiMimi.room_relabel_parameter`** &mdash; *Method*.



Connect a gradient to another component: change the variable component and name to another component.

<a id='OptiMimi.roomdiagonal-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function}' href='#OptiMimi.roomdiagonal-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function}'>#</a>
**`OptiMimi.roomdiagonal`** &mdash; *Method*.



```
roomdiagonal(model, component, variable, parameter, gen)
```

Fill in just the diagonal, assuming $\frac{dv_i}{dp_j} = 0$, if $i <> j$.  Requires that the dimensions of `variable` and `parameter` are the same.  `gen` called for each combination of indices.

$\begin{array}{ccc} p_1 & p_2 & \cdots \end{array}$
$\begin{array}{c}
    v_1 \\
    v_2 \\
    \vdots
    \end{array}\left(\begin{array}{ccc}
        g(1) & 0 & \cdots \\
        0 & g(2) & \cdots \\
        \vdots & \vdots & \ddots
        \end{array}\right)$
**Arguments**

  * `model::Model`: The model containing `component`.
  * `component::Symbol`: The component containing `variable` and `parameter`.
  * `variable::Symbol`: The outcome variable, corresponding to matrix rows.
  * `parameter::Symbol`: The optimization parameter, corresponding to matrix columns.
  * `gen::Function`: The function generating gradient values.

<a id='OptiMimi.roomintersect-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function}' href='#OptiMimi.roomintersect-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function}'>#</a>
**`OptiMimi.roomintersect`** &mdash; *Method*.



```
roomintersect(model, component, variable1, variable2, gen)
```

Generate a room at the intersection of two variables.

Call `gen` for each shared index, passing an array to be filled for unshared indexes. This version assumes both variables come from the same component.

See `matrixintersect` for the matrix generation logic.

<a id='OptiMimi.roomintersect-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Symbol,Function}' href='#OptiMimi.roomintersect-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Symbol,Function}'>#</a>
**`OptiMimi.roomintersect`** &mdash; *Method*.



```
roomintersect(model, component, variable1, component2, variable2, gen)
```

Generate a room at the intersection of two variables.

Call `gen` for each shared index, passing an array to be filled for unshared indexes. This version allows the variables to come from different component.

See `matrixintersect` for the matrix generation logic.

<a id='OptiMimi.roomsingle-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function}' href='#OptiMimi.roomsingle-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function}'>#</a>
**`OptiMimi.roomsingle`** &mdash; *Method*.



Fill in every element

<a id='OptiMimi.setconstraintoffset!-Tuple{OptiMimi.LinearProgrammingHouse,OptiMimi.LinearProgrammingHall}' href='#OptiMimi.setconstraintoffset!-Tuple{OptiMimi.LinearProgrammingHouse,OptiMimi.LinearProgrammingHall}'>#</a>
**`OptiMimi.setconstraintoffset!`** &mdash; *Method*.



Set offset values from a `LinearProgrammingHall`.  See the other `setconstraintoffset!`

<a id='OptiMimi.setconstraintoffset!-Tuple{OptiMimi.LinearProgrammingHouse,Symbol,Symbol,Array{Float64,1}}' href='#OptiMimi.setconstraintoffset!-Tuple{OptiMimi.LinearProgrammingHouse,Symbol,Symbol,Array{Float64,1}}'>#</a>
**`OptiMimi.setconstraintoffset!`** &mdash; *Method*.



```
setconstraintoffset!(house, component, variable, f)
```

Set offset values within a `LinearProgrammingHouse`.

**Arguments**

  * `house::LinearProgrammingHouse`: The house to set values within.
  * `component::Symbol`: The component for the constraint variable.
  * `variable::Symbol`: The variable for the constraint.
  * `f::Vector{Float64}`: The values, with all dimensions collapsed into a single vector.

<a id='OptiMimi.varlengths' href='#OptiMimi.varlengths'>#</a>
**`OptiMimi.varlengths`** &mdash; *Function*.



Return the total span occupied by each variable or parameter.

<a id='OptiMimi.findinfeasiblepairhelper-Tuple{OptiMimi.LinearProgrammingHouse,Any,Any,Any,Any}' href='#OptiMimi.findinfeasiblepairhelper-Tuple{OptiMimi.LinearProgrammingHouse,Any,Any,Any,Any}'>#</a>
**`OptiMimi.findinfeasiblepairhelper`** &mdash; *Method*.



Look for the row that makes the problem feasible

<a id='OptiMimi.getdimnames-Tuple{Mimi.Model,Symbol,Symbol}' href='#OptiMimi.getdimnames-Tuple{Mimi.Model,Symbol,Symbol}'>#</a>
**`OptiMimi.getdimnames`** &mdash; *Method*.



Return the symbols representing each of the dimensions for this variable or parameter.

<a id='OptiMimi.getdims-Tuple{Mimi.Model,Symbol,Symbol}' href='#OptiMimi.getdims-Tuple{Mimi.Model,Symbol,Symbol}'>#</a>
**`OptiMimi.getdims`** &mdash; *Method*.



Return a vector of the indices defining the parameter or variable.

<a id='OptiMimi.matrixdiagonal-Tuple{Array{Int64,1},Function}' href='#OptiMimi.matrixdiagonal-Tuple{Array{Int64,1},Function}'>#</a>
**`OptiMimi.matrixdiagonal`** &mdash; *Method*.



```
matrixdiagonal(dims, gen)
```

Creates a matrix of dimensions $\prod \text{dims}_i$.  Call the generate function, `gen` for all indices along the diagonal.  All combinations of indices will be called, since the "diagonal" part is between the rows and the columns.

**Arguments**

  * `dims::Vector{Int64}`: The multiple dimensions collapsed into both the rows and columns.
  * `gen::Function`: A function called with an argument for each dimension.

**Examples**

```julia
using OptiMimi

dims = [3, 2]
A = OptiMimi.matrixdiagonal(dims, (ii, jj) -> 1)
sum(A)

# output
6.0
```

<a id='OptiMimi.matrixintersect-Tuple{Array{Int64,1},Array{Int64,1},Array{Symbol,1},Array{Symbol,1},Function}' href='#OptiMimi.matrixintersect-Tuple{Array{Int64,1},Array{Int64,1},Array{Symbol,1},Array{Symbol,1},Function}'>#</a>
**`OptiMimi.matrixintersect`** &mdash; *Method*.



```
matrixintersect(rowdims, coldims, rowdimnames, coldimnames, gen)
```

Call the `gen` function with all combinations of the shared indices. Shared dimensions must come in the same order and at the end of the dimensions lists.

**Arguments**

  * `rowdims::Vector{Int64}` and `coldims::Vector{Int64}`: the number of elements in each dimension of the row and column variables.
  * `rowdimnames::Vector{Symbol}` and `coldimnames::Vector{Symbol}`: The names of the dimensions; note that the last 1 or more dimensions should be the same in both of these lists.
  * `gen::Function`: a function of arguments (`A`, `indices`...), where `A` is a matrix of the un-shared indices and `indices` are multiple indexes describing all shared indices.

**Examples**

```julia
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

<a id='OptiMimi.matrixsingle-Tuple{Array{Int64,1},Array{Int64,1},Any}' href='#OptiMimi.matrixsingle-Tuple{Array{Int64,1},Array{Int64,1},Any}'>#</a>
**`OptiMimi.matrixsingle`** &mdash; *Method*.



Call the generate function for every element.

<a id='OptiMimi.toindex-Tuple{Int64,Array{Int64,1}}' href='#OptiMimi.toindex-Tuple{Int64,Array{Int64,1}}'>#</a>
**`OptiMimi.toindex`** &mdash; *Method*.



Translate an offset value (+1) to an index vector.

<a id='OptiMimi.vectorsingle-Tuple{Array{Int64,1},Any}' href='#OptiMimi.vectorsingle-Tuple{Array{Int64,1},Any}'>#</a>
**`OptiMimi.vectorsingle`** &mdash; *Method*.



Construct a matrix with the given dimensions, calling gen for each element.

