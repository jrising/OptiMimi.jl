module OptiMimi

using NLopt
using ForwardDiff, DiffBase
using MathProgBase
using Compat

import Mimi: Model

export problem, solution, unaryobjective, objevals, setparameters, nameindexes, sensitivity, uncertainproblem

include("registerdiff.jl")
include("matrixconstraints.jl")
include("linproghouse.jl")
include("makeroom.jl")
include("metamimi.jl")
include("uncertainty.jl")

allverbose = false
objevals = 0

type OptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    opt::Opt
    constraints::Vector{Function}
end

type BlackBoxOptimizationProblem{T}
    algorithm::Symbol
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    objective::Function
    lowers::Vector{T}
    uppers::Vector{T}
end

type LinprogOptimizationProblem{T}
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    objective::Function
    objectiveconstraints::Vector{Function}
    matrixconstraints::Vector{MatrixConstraintSet}
    exlowers::Vector{T}
    exuppers::Vector{T}
end

BlackBoxAlgorithms = [:separable_nes, :xnes, :dxnes]

"""Returns (ii, len, isscalar) with the index of each symbol and its length."""
function nameindexes(model::Model, components::Vector{Symbol}, names::Vector{Symbol})
    for ii in 1:length(components)
        dims = getdims(model, components[ii], names[ii])
        if length(dims) == 0
            produce((ii, 1, true)) # It's a scalar
        else
            produce((ii, prod(dims), false))
        end
    end
end

"""Set parameters in a model."""
function setparameters(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, xx::Vector)
    startindex = 1
    for (ii, len, isscalar) in @task nameindexes(model, components, names)
        if isscalar
            setfield!(get(model.mi).components[components[ii]].Parameters, names[ii], xx[startindex])
        else
            shape = getdims(model, components[ii], names[ii])
            reshaped = reshape(collect(model.numberType, xx[startindex:(startindex+len - 1)]), tuple(shape...))
            setfield!(get(model.mi).components[components[ii]].Parameters, names[ii], reshaped)
        end
        startindex += len
    end
end

function sensitivity(model::Model, component::Symbol, parameter::Symbol, objective::Function, points::Vector{Float64})
    results = []
    for point in points
        setparameters(model, [component], [parameter], [point])
        run(model)
        push!(results, objective(model))
    end

    results
end

"""Generate the form of objective function used by the optimization, taking parameters rather than a model."""
function unaryobjective(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, objective::Function)
    function my_objective(xx::Vector)
        if allverbose
            println(xx)
        end

        global objevals
        objevals += 1

        setparameters(model, components, names, xx)
        run(model)
        objective(model)
    end

    my_objective
end

"""Create an NLopt-style objective function which does not use its grad argument."""
function gradfreeobjective(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, objective::Function)
    myunaryobjective = unaryobjective(model, components, names, objective)
    function myobjective(xx::Vector, grad::Vector)
        myunaryobjective(xx)
    end

    myobjective
end

"""Create an NLopt-style objective function which computes an autodiff gradient."""
function autodiffobjective(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, objective::Function)
    myunaryobjective = unaryobjective(model, components, names, objective)
    function myobjective(xx::Vector, grad::Vector)
        out = DiffBase.GradientResult(xx)
        ForwardDiff.gradient!(out, myunaryobjective, xx)
        if any(isnan(DiffBase.gradient(out)))
            error("objective gradient is NaN")
        end
        copy!(grad, DiffBase.gradient(out))
        DiffBase.value(out)
    end

    myobjective
end

"""Create a 0 point."""
function make0(model::Model, components::Vector{Symbol}, names::Vector{Symbol})
    initial = Float64[]
    for (ii, len, isscalar) in @task nameindexes(model, components, names)
        append!(initial, [0. for jj in 1:len])
    end

    initial
end

"""Expand parameter constraints to full vectors for every numerical parameter."""
function expandlimits{T<:Real}(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{T}, uppers::Vector{T})
    my_lowers = T[]
    my_uppers = T[]

    ## Replace with eachname
    totalvars = 0
    for (ii, len, isscalar) in @task nameindexes(model, components, names)
        append!(my_lowers, [lowers[ii] for jj in 1:len])
        append!(my_uppers, [uppers[ii] for jj in 1:len])
        totalvars += len
    end

    my_lowers, my_uppers, totalvars
end

"""Setup an optimization problem."""
function problem{T<:Real}(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{T}, uppers::Vector{T}, objective::Function; constraints::Vector{Function}=Function[], algorithm::Symbol=:LN_COBYLA_OR_LD_MMA)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)

    if algorithm == :GUROBI_LINPROG
        # Make no changes to objective!
    elseif model.numberType == Number
        if algorithm == :LN_COBYLA_OR_LD_MMA
            algorithm = :LD_MMA
        end
        if string(algorithm)[2] == 'N'
            warn("Model is autodifferentiable, but optimizing using a derivative-free algorithm.")
            myobjective = gradfreeobjective(model, components, names, objective)
        else
            println("Using AutoDiff objective.")
            myobjective = autodiffobjective(model, components, names, objective)
        end
    else
        if algorithm == :LN_COBYLA_OR_LD_MMA
            algorithm = :LN_COBYLA
        elseif string(algorithm)[2] == 'D'
            warn("Model is non-differentiable, but requested a gradient algorithm; instead using LN_COBYLA.")
            algorithm = :LN_COBYLA
        end

        myobjective = gradfreeobjective(model, components, names, objective)
    end

    if algorithm == :GUROBI_LINPROG
        LinprogOptimizationProblem(model, components, names, objective, constraints, MatrixConstraintSet[], my_lowers, my_uppers)
    else
        if algorithm in BlackBoxAlgorithms
            if length(constraints) > 0
                warn("Functional constraints not supported for BBO algorithms.")
            end

            BlackBoxOptimizationProblem(algorithm, model, components, names, myobjective, my_lowers, my_uppers)
        else
            opt = Opt(algorithm, totalvars)
            lower_bounds!(opt, my_lowers)
            upper_bounds!(opt, my_uppers)
            xtol_rel!(opt, 1e-6)

            max_objective!(opt, myobjective)

            for constraint in constraints
                let this_constraint = constraint
                    function my_constraint(xx::Vector, grad::Vector)
                        setparameters(model, components, names, xx)
                        this_constraint(model)
                    end

                    inequality_constraint!(opt, my_constraint)
                end
            end

            OptimizationProblem(model, components, names, opt, constraints)
        end
    end
end

"""Solve an optimization problem."""
function solution(optprob::OptimizationProblem, generator::Function; maxiter=Inf, verbose=false)
    global allverbose
    allverbose = verbose

    if verbose
        println("Selecting an initial point.")
    end

    attempts = 0
    initial = []
    valid = false
    while attempts < maxiter
        initial = generator()

        setparameters(optprob.model, optprob.components, optprob.names, initial)

        valid = true
        for constraint in optprob.constraints
            if constraint(optprob.model) >= 0
                valid = false
                break
            end
        end

        if valid
            break
        end

        attempts += 1
        if attempts % 1000 == 0
            println("Could not find initial point after $attempts attempts.")
        end
    end

    if !valid
        throw(DomainError("Could not find a valid initial value."))
    end

    if verbose
        println("Optimizing...")
    end
    if maxiter < Inf
        maxeval!(optprob.opt, maxiter)
    end
    (minf,minx,ret) = optimize(optprob.opt, initial)

    (minf, minx)
end

"""Solve an optimization problem."""
function solution(optprob::BlackBoxOptimizationProblem; maxiter=Inf, verbose=false)
    global allverbose
    allverbose = verbose

    if verbose
        println("Optimizing...")
    end
    res = bboptimize(p -> -optprob.objective(p, 0 .* p); SearchRange=collect(zip(optprob.lowers, optprob.uppers)), Method=optprob.algorithm)

    (best_fitness(res), best_candidate(res))
end


"""Setup an optimization problem."""
function problem{T<:Real}(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{T}, uppers::Vector{T}, objective::Function, objectiveconstraints::Vector{Function}, matrixconstraints::Vector{MatrixConstraintSet})
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)
    LinprogOptimizationProblem(model, components, names, objective, objectiveconstraints, matrixconstraints, my_lowers, my_uppers)
end

"""Solve an optimization problem."""
function solution(optprob::LinprogOptimizationProblem, verbose=false)
    global allverbose
    allverbose = verbose

    initial = make0(optprob.model, optprob.components, optprob.names)

    if verbose
        println("Optimizing...")
    end

    if optprob.model.numberType == Number
        myobjective = unaryobjective(optprob.model, optprob.components, optprob.names, optprob.objective)
        f, b, A = lpconstraints(optprob.model, optprob.components, optprob.names, myobjective, objectiveconstraints)
    else
        f, b, A = lpconstraints(optprob.model, optprob.components, optprob.names, optprob.objective, optprob.objectiveconstraints)
    end

    f, b, A = combineconstraints(f, b, A, optprob.model, optprob.components, optprob.names, optprob.matrixconstraints)
    exlowers, exuppers = combinelimits(optprob.exlowers, optprob.exuppers, optprob.model, optprob.components, optprob.names, optprob.matrixconstraints)

    # Use -f, because linprog *minimizes* objective
    @time sol = linprog(-f, A, '<', b, optprob.exlowers, optprob.exuppers)

    sol.sol
end

end # module
