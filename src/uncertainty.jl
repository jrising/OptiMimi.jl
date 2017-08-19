# We use genetic algorithms (https://github.com/WestleyArgentum/GeneticAlgorithms.jl) because:
# 1. The size of the state space >> parameter space.  So just explore parameter space.
# 2. We cannot run from the last period reliably without knowing the state space.

using BlackBoxOptim

type UncertainOptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    objective::Function
    lowers::Vector{Float64}
    uppers::Vector{Float64}
end

"""Setup an optimization over Monte Carlo uncertainty."""
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function, mcperlife::Int64)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)

    my_unaryobjective = unaryobjective(model, components, names, objective)
    function my_objective(parameters::Vector{Float64})
        total = 0
        for iter in 1:mcperlife
            montecarlo(model)
            total += my_unaryobjective(parameters)
        end

        -total / mcperlife
    end

    UncertainOptimizationProblem(model, components, names, my_objective, my_lowers, my_uppers)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem)
    res = bboptimize(optprob.objective; SearchRange=collect(zip(optprob.lowers, optprob.uppers)), MaxFuncEvals=10000 * 5)
    best_candidate(res)
end

