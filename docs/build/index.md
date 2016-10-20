
<a id='OptiMimi-Documentation-1'></a>

# OptiMimi Documentation


OptiMimi provides a simplified interface for finding optimal parameter values for Mimi models.  The package contains two major interfaces:


  * An abstraction around general nonlinear parameter optimization, using simulation runs.
  * A collection of helper tools that work with Mimi models to generate linear programming matrices.


<a id='Installation-1'></a>

## Installation


At the Julia REPL:


```julia
    Pkg.add("OptiMimi")
    Pkg.checkout("OptiMimi")
```


The second line ensures that you have the most up-to-date version.


<a id='Usage-1'></a>

## Usage


For the nonlinear parameter optimization system, see [the ReadMe](https://github.com/jrising/OptiMimi.jl/blob/master/README.md).


For the linear programming system, see the [Linear Programming Guide](linproghouse.md#Linear-Programming-Guide-1).


<a id='Index-1'></a>

## Index

- [`OptiMimi.LinearProgrammingHall`](linproghouse.md#OptiMimi.LinearProgrammingHall)
- [`OptiMimi.LinearProgrammingHouse`](linproghouse.md#OptiMimi.LinearProgrammingHouse)
- [`OptiMimi.LinearProgrammingRoom`](linproghouse.md#OptiMimi.LinearProgrammingRoom)
- [`OptiMimi.LinearProgrammingShaft`](linproghouse.md#OptiMimi.LinearProgrammingShaft)
- [`OptiMimi.findinfeasiblepair`](linproghouse.md#OptiMimi.findinfeasiblepair-Tuple{OptiMimi.LinearProgrammingHouse,Any})
- [`OptiMimi.findinfeasiblepairhelper`](linproghouse.md#OptiMimi.findinfeasiblepairhelper-Tuple{OptiMimi.LinearProgrammingHouse,Any,Any,Any,Any})
- [`OptiMimi.fromindex`](linproghouse.md#OptiMimi.fromindex-Tuple{Array{Int64,1},Array{Int64,1}})
- [`OptiMimi.getconstraintoffset`](linproghouse.md#OptiMimi.getconstraintoffset-Tuple{OptiMimi.LinearProgrammingHouse,Symbol,Symbol})
- [`OptiMimi.getconstraintsolution`](linproghouse.md#OptiMimi.getconstraintsolution-Tuple{Any,Any,Any})
- [`OptiMimi.getdimnames`](linproghouse.md#OptiMimi.getdimnames-Tuple{Mimi.Model,Symbol,Symbol})
- [`OptiMimi.getdims`](linproghouse.md#OptiMimi.getdims-Tuple{Mimi.Model,Symbol,Symbol})
- [`OptiMimi.getparametersolution`](linproghouse.md#OptiMimi.getparametersolution-Tuple{OptiMimi.LinearProgrammingHouse,Array{Float64,1},Symbol})
- [`OptiMimi.hall_relabel`](linproghouse.md#OptiMimi.hall_relabel-Tuple{OptiMimi.LinearProgrammingHall,Symbol,Symbol,Symbol})
- [`OptiMimi.matrixdiagonal`](linproghouse.md#OptiMimi.matrixdiagonal-Tuple{Array{Int64,1},Function})
- [`OptiMimi.matrixintersect`](linproghouse.md#OptiMimi.matrixintersect-Tuple{Array{Int64,1},Array{Int64,1},Array{Symbol,1},Array{Symbol,1},Function})
- [`OptiMimi.matrixsingle`](linproghouse.md#OptiMimi.matrixsingle-Tuple{Array{Int64,1},Array{Int64,1},Any})
- [`OptiMimi.room_relabel`](linproghouse.md#OptiMimi.room_relabel-Tuple{OptiMimi.LinearProgrammingRoom,Symbol,Symbol,Symbol})
- [`OptiMimi.room_relabel_parameter`](linproghouse.md#OptiMimi.room_relabel_parameter-Tuple{OptiMimi.LinearProgrammingRoom,Symbol,Symbol,Symbol})
- [`OptiMimi.roomdiagonal`](linproghouse.md#OptiMimi.roomdiagonal-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function})
- [`OptiMimi.roomintersect`](linproghouse.md#OptiMimi.roomintersect-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Symbol,Function})
- [`OptiMimi.roomintersect`](linproghouse.md#OptiMimi.roomintersect-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function})
- [`OptiMimi.roomsingle`](linproghouse.md#OptiMimi.roomsingle-Tuple{Mimi.Model,Symbol,Symbol,Symbol,Function})
- [`OptiMimi.setconstraintoffset!`](linproghouse.md#OptiMimi.setconstraintoffset!-Tuple{OptiMimi.LinearProgrammingHouse,Symbol,Symbol,Array{Float64,1}})
- [`OptiMimi.setconstraintoffset!`](linproghouse.md#OptiMimi.setconstraintoffset!-Tuple{OptiMimi.LinearProgrammingHouse,OptiMimi.LinearProgrammingHall})
- [`OptiMimi.toindex`](linproghouse.md#OptiMimi.toindex-Tuple{Int64,Array{Int64,1}})
- [`OptiMimi.varlengths`](linproghouse.md#OptiMimi.varlengths)
- [`OptiMimi.vectorsingle`](linproghouse.md#OptiMimi.vectorsingle-Tuple{Array{Int64,1},Any})

