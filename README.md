# OptiMimi: Optimization of Mimi models

OptiMimi provides a simplified interface for finding optimal parameter
values for Mimi models (https://github.com/anthofflab/Mimi.jl).  The
core interface consists of `problem` to define the optimization
problem, and calls to `solution` to solve it.

The package provides two major approaches for performing optimization:
General and Linear programming.  The General approach allows takes
full models and applies non-linear optimization techniques to optimize
parameters within them.  The Linear programming approach allows models
to define linear operation matrices which represent the computations
they perform in Mimi's `run_timestep` function.  In addition, OptiMimi
offers a way to automatically generate these matrices.

OptiMimi supports autodifferentiation using ForwardDiff.  To use it,
the Model must be created with the optional `autodiffable` set to
`true`, and all components must be created using the `@defcompo`
macro, instead of `@defcomp`.  If these are not used, OptiMimi will
fall back on derivative-free algorithms.

The General approach in OptiMimi can use algorithms in both NLopt and
BlackBoxOptim.  The Linear programming approach uses any solver
supported by MathProgBase.

## Constructing an optimization problem (General approach)

Setup an optimization problem with the `problem` function:
```
problem(model, names, lowers, uppers, objective; constraints, algorithm)
```

* `model` is a Mimi model, with some parameters intended for optimization.
* `components` and `names` are lists of the parameter names to be optimized, and must be the same length.
* `lowers` and `uppers` are a list of lower and upper bounds and must be the same length as `names`; the same bounds are used for all values within a parameter array.
* `objective` is a function that takes a Mimi Model, with parameters set and fully executed, and returns a value to be maximized.
* `constraints` (optional) is a vector of inequality constraint functions; each takes a Mimi Model, with parameters set (but not necessarily executed), and should return < 0 when the constraint is satisfied.
* `algorithm` (optional) is a symbol, currently chosen from the NLopt algorithms.

The return value is an object of the OptimizationProblem type, to be passed to `solution`.

### Example:

Start by creating a Mimi model and ensuring that it runs with all
parameters set.  In the example below, `my_model` is a model with an
agriculture component, in which N regions are evaluated in a single
timestep to consume energy and produce corn.

The optimization maximizes economic value, trading off the value of
the corn against the cost of the energy for fertilizer.  We also add a
constraint that the total fertilizer cannot be more than 1 million kg,
to reduce environmental impacts.

```
using OptiMimi

# Prices of goods
p_F = 0.25  # the global price of food (per kg of corn)
p_E = 0.4   # the global price of fuel (per kWh)

# Objective to maximize economic output
function objective(model::Model)
    sum(my_model[:agriculture, :cornproduction] * p_F - my_model[:agriculture, :cornenergyuse] * p_E)
end

constraints = [model -> sum(model.components[:agriculture].Parameters.fertilizer) - 1e6]

# Setup the optimization
optprob = problem(my_model, [:agriculture], [:fertilizer], [0.], [1e6], objective, constraints=constraints)
```

Note that (1) the objective function is provided with the prepared
model, not with the raw initialization values, and (2) even though
there are N values to be set and optimized over in the `fertilizer`
parameter, the lower and upper bounds are only specified once.

## Solving the optimization problem

The optimization problem, returned by `problem` is solved by `solution`:
```
solution(optprob, generator; maxiter, verbose)
```

* `optprob` is the result of the `problem` function.
* `generator` is a function of no arguments, which returns a full set of parameter values, with values concatenated across parameters in the order of `names` above.  This should generally be stochastic, and if the specified model fails the constraints then `generator` will be called again until it succeeds.
* `maxiter` (optional) is the maximum number of iterations for the optimization; currently it only is used for the maximum number of times that `generator` will be called.
* `verbose` (optional) is a boolean to specify if status messages should be printed.

The return value is a tuple of the maximum found objective value, and
the concatenated collection of model parameters that produced it.

### Example:

Continuing the example above, we solve the optimization problem:

```
(maxf, maxx) = solution(optprob, () -> [0. for i in 1:5])

println(maxf)
println(maxx)
```

Our generator function can only generate a single initial condition: all 0's.

## Constructing an optimization problem (Linear programming approach)

Linear programming allows for vastly faster optimizations, so long as
the constraints and objective can be translated into linear algebra
relationships.  See https://en.wikipedia.org/wiki/Linear_programming
for details.

Within OptiMimi, the large matrix and vectors which define the linear
programming constraints and objective are developed in segments.
In line with the model organization structure of Mimi, the
computations which relate variables to parameters are kept separate,
and organized with each component.  These computations can be used
directly, if variables of a component are constrained and the
parameters of that same component are optimized over; or, they can be
connected across multiple components.  Rather than defining the entire
matrix at once, OptiMimi allows segments of the matrix specific to
each component to be defined separately.

Segments of the constraint matrix are encapsulated in
`LinearProgrammingRoom` objects, which combine a sparse matrix with
information about a single model parameter and variable.  A column
vector, combined with information about a single model parameter or
variable, is encapsulated in a `LinearProgrammingHall` or
`LinearProgrammingShaft` (its transpose) object.  A
`LinearProgrammingHouse` contains the set of all matrices for the
optimization.

There are a variety of functions available which create these objects
or manipulate them.  Some of the most commonly used operations are
below:

 * `roomdiagonal`: creates a room for a variable which is a direct
   scaling of a parameter.
 * `*`: Allows two rooms to be multiplied together, which corresponds
   to "connecting" the variable of the first as the parameter of the
   second; or allows a room and a hall or shaft to be multiplied so
   which results in a weighted sum of the variables in the room (e.g.,
   as an objective).
 * `room_relabel`: Relabel the variable of a room so that it
   corresponds to the parameter name of another room which it is to
   connected to (multiplied with).

The usual process for setting up a linear programming problem is as
follows:

1. Mimi components are written as normal, with `@defcomp` calls, and
   the Mimi model is constructed and external parameters are
   initialized.
   
2. Individual functions are specified for each component describing
   the gradient of a variable with respect to a parameter, using the
   `room` functions or the automated option in `makeroom.jl`.  The
   naming typically is as follows:
   ```grad_COMPONENT_VARIABLE_PARAMETER(model)```
   ```constraintoffset_COMPONENT_VARIABLE(model)```
   
2. A `LinearProgrammingHosue` is constructed, specifying the
   optimization parameters and constraint variables, like so:
   ```house = LinearProgrammingHouse(model, paramcomponents,
   parameters, constcomponents, variables)```
   
3. The objective is specified with a `setobjective!` call, often as a
   sum over variable specified by a gradient function, e.g.,
   ```setobjective!(house,
   -varsum(grad_COMPONENT_cost_PARAMETER(model)))```
   
4. Constraints are specified with `setconstraint!` and
   `setconstraintoffset!` calls.  In all cases, the relationship must
   be specified so that `variable < offset`.  This looks like,
   ```setconstraint!(house,
   grad_COMPONENT_VARIABLE_PARAMETER(model))```
   ```setconstraintoffset!(house,
   constraintoffset_COMPONENT_VARIABLE(model))```
   
5. The optimization is performed, using any solver supported by
   `MathProgBase` and the `houseoptimize` function.  For example:
```
using MathProgBase
using Gurobi
solver = GurobiSolver()

sol = houseoptimize(house, solver)
```

6. The result is studied using the `summarizeparameters`, or, if the
   optimzation failed, `findinfeasiblepair`.
