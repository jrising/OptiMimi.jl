# OptiMimi Documentation

OptiMimi provides a simplified interface for finding optimal parameter
values for Mimi models.  The package contains two major interfaces:

* An abstraction around general nonlinear parameter optimization, using simulation runs.
* A collection of helper tools that work with Mimi models to generate linear programming matrices.

## Installation

At the Julia REPL:

```julia
    Pkg.add("OptiMimi")
    Pkg.checkout("OptiMimi")
```

The second line ensures that you have the most up-to-date version.

## Usage

For the nonlinear parameter optimization system, see [the ReadMe](https://github.com/jrising/OptiMimi.jl/blob/master/README.md).

For the linear programming system, see the [Linear Programming Guide](@ref).

## Index

```@index
Pages = ["linproghouse.md"]
```

