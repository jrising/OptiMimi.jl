<a id='OptiMimi.nameindexes-Tuple{Mimi.Model,Array{Symbol,1}}' href='#OptiMimi.nameindexes-Tuple{Mimi.Model,Array{Symbol,1}}'>#</a>
**`OptiMimi.nameindexes`** &mdash; *Method*.



Returns (ii, len, isscalar) with the index of each symbol and its length.

<a id='OptiMimi.problem-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Array{T<:Real,1},Array{T<:Real,1},Function,Array{Function,1},Array{OptiMimi.MatrixConstraintSet,1}}' href='#OptiMimi.problem-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Array{T<:Real,1},Array{T<:Real,1},Function,Array{Function,1},Array{OptiMimi.MatrixConstraintSet,1}}'>#</a>
**`OptiMimi.problem`** &mdash; *Method*.



Setup an optimization problem.

<a id='OptiMimi.problem-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Array{T<:Real,1},Array{T<:Real,1},Function}' href='#OptiMimi.problem-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Array{T<:Real,1},Array{T<:Real,1},Function}'>#</a>
**`OptiMimi.problem`** &mdash; *Method*.



Setup an optimization problem.

<a id='OptiMimi.setparameters-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Array{T,1}}' href='#OptiMimi.setparameters-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Array{T,1}}'>#</a>
**`OptiMimi.setparameters`** &mdash; *Method*.



Set parameters in a model.

<a id='OptiMimi.solution' href='#OptiMimi.solution'>#</a>
**`OptiMimi.solution`** &mdash; *Function*.



Solve an optimization problem.

<a id='OptiMimi.solution-Tuple{OptiMimi.OptimizationProblem,Function}' href='#OptiMimi.solution-Tuple{OptiMimi.OptimizationProblem,Function}'>#</a>
**`OptiMimi.solution`** &mdash; *Method*.



Solve an optimization problem.

<a id='OptiMimi.unaryobjective-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Function}' href='#OptiMimi.unaryobjective-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Function}'>#</a>
**`OptiMimi.unaryobjective`** &mdash; *Method*.



Generate the form of objective function used by the optimization, taking parameters rather than a model.

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

$egin{array}{ccc} p_1 & p_2 & cdots nd{array}$
$egin{array}{c}
    v_1 \
    v_2 \
    dots
    nd{array}left(egin{array}{ccc}
        g(1) & 0 & cdots \
        0 & g(2) & cdots \
        dots & dots & ddots
        nd{array}ight)$
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

<a id='OptiMimi.autodiffobjective-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Function}' href='#OptiMimi.autodiffobjective-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Function}'>#</a>
**`OptiMimi.autodiffobjective`** &mdash; *Method*.



Create an NLopt-style objective function which computes an autodiff gradient.

<a id='OptiMimi.gradfreeobjective-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Function}' href='#OptiMimi.gradfreeobjective-Tuple{Mimi.Model,Array{Symbol,1},Array{Symbol,1},Function}'>#</a>
**`OptiMimi.gradfreeobjective`** &mdash; *Method*.



Create an NLopt-style objective function which does not use its grad argument.

<a id='OptiMimi.make0-Tuple{Mimi.Model,Array{Symbol,1}}' href='#OptiMimi.make0-Tuple{Mimi.Model,Array{Symbol,1}}'>#</a>
**`OptiMimi.make0`** &mdash; *Method*.



Create a 0 point.

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

Creates a matrix of dimensions $prod 	ext{dims}_i$.  Call the generate function, `gen` for all indices along the diagonal.  All combinations of indices will be called, since the "diagonal" part is between the rows and the columns.

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

