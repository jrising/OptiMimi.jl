# We use genetic algorithms (https://github.com/WestleyArgentum/GeneticAlgorithms.jl) because:
# 1. The size of the state space >> parameter space.  So just explore parameter space.
# 2. We cannot run from the last period reliably without knowing the state space.

using GeneticAlgorithms

module MimiGA

using GeneticAlgorithms
using Distributions

paramsize = 0
lowers = Array(Float64, 0)
uppers = Array(Float64, 0)
objective = (parameters) -> 0
iterations = 0

function prepare(lows::Vector{Float64}, his::Vector{Float64}, obj::Function)
    global paramsize, lowers, uppers, objective, iterations

    paramsize = length(lows)
    lowers = lows
    uppers = his
    objective = obj
    iterations = 0
end

type ParameterSet <: Entity
    parameters::Vector{Float64}
    fitness

    ParameterSet() = new(Array(Float64, paramsize), nothing)
    ParameterSet(parameters::Vector{Float64}) = new(parameters, nothing)
end

function create_entity(num)
    # for simplicity sake, we will use uniform values between the bounds
    ParameterSet(rand(paramsize) .* (uppers - lowers) + lowers)
end

function fitness(ent)
    obj = objective(ent.parameters)
    if ent.fitness != nothing
        ent.fitness * .9 + obj * .1
    else
        obj * .1
    end
end

function group_entities(pop)
    global iterations

    iterations = iterations + 1
    if iterations % 100 == 0
        println(pop[1])
	if iterations >= 10000
	    return
	end
    end

    # simple naive groupings that pair the best entitiy with every other
    dist = TriangularDist(1, length(pop)+1, 1)
    for i in 1:length(pop)
        produce([floor(Int, rand(dist)), floor(Int, rand(dist))])
    end
end

function crossover(group)
    child = ParameterSet()

    # grab each element from a random parent
    num_parents = length(group)
    for i in 1:length(group[1].parameters)
        parent = (rand(UInt) % num_parents) + 1
        child.parameters[i] = group[parent].parameters[i]
    end

    child.fitness = (group[1].fitness + group[2].fitness) / 2

    child
end

function mutate(ent)
    rand(Float64) < 0.8 && return

    randii = rand(UInt) % paramsize + 1
    current = (ent.parameters[randii] - lowers[randii]) / (uppers[randii] - lowers[randii])
    dist = Beta(current * iterations + 1, (1 - current) * iterations + 1)
    ent.parameters[randii] = rand(dist) * (uppers[randii] - lowers[randii]) + lowers[randii]
end

end

type UncertainOptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    objective::Function
    lowers::Vector{Float64}
    uppers::Vector{Float64}
end


"""Setup an optimization over Monte Carlo uncertainty."""
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)

    my_unaryobjective = unaryobjective(model, components, names, objective)
    function my_objective(parameters::Vector{Float64})
        montecarlo(model)
        my_unaryobjective(parameters)
    end

    UncertainOptimizationProblem(model, components, names, my_objective, my_lowers, my_uppers)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem)
    MimiGA.prepare(optprob.lowers, optprob.uppers, optprob.objective)

    gamodel = runga(MimiGA; initial_pop_size = 16)

    GeneticAlgorithms.population(gamodel)  # the the latest population when the GA exited
end

