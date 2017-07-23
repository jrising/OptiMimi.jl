"""
We use genetic algorithms (https://github.com/WestleyArgentum/GeneticAlgorithms.jl) because:

1. The size of the state space >> parameter space.  So just explore parameter space.
2. We cannot run from the last period reliably without knowing the state space.
"""
using GeneticAlgorithms

module MimiGA

paramsize = 0
lowers = Array(Float64, 0)
uppers = Array(Float64, 0)
objective(parameters) = 0

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
    objective(ent.parameters)
end

function group_entities(pop)
    # simple naive groupings that pair the best entitiy with every other
    for i in 1:length(pop)
        produce([1, i])
    end
end

function crossover(group)
    child = ParameterSet()

    # grab each element from a random parent
    num_parents = length(group)
    for i in 1:length(group[1].parameters)
        parent = (rand(Uint) % num_parents) + 1
        child.parameters[i] = group[parent].parameters[i]
    end

    child
end

function mutate(ent)
    rand_element = rand(Uint) % paramsize + 1
    ent.parameters[rand_element] = rand() .* (uppers - lowers) + lowers
end

end

type UncertainOptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    objective::Function
    lowers::Vector{Float64}
    uppers::Vector{Float64}
    montecarlo::Function
end


"""Setup an optimization over Monte Carlo uncertainty."""
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)

    my_unaryobjective = unaryobjective(model, components, names, objective)
    function my_objective(parameters::Vector{Float64})
        montecarlo(model)
        my_unaryobjective(parameters)
    end

    UncertainOptimizationProblem(model, components, names, my_objective, my_lower, my_uppers)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem, generator::Function; maxiter=Inf, verbose=false)
    MimiGA.paramsize = length(optprob.lowers)
    MimiGA.lowers = optprob.lowers
    MimiGA.uppers = optprob.uppers
    MimiGA.objective = optprob.objective

    gamodel = runga(MimiGA; initial_pop_size = 16)

    population(gamodel)  # the the latest population when the GA exited
end

