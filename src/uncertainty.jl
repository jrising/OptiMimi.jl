## We use Monte Carlo methods, because: 1. The size of the state space
## >> parameter space.  So just explore parameter space.  2. We cannot
## run from the last period reliably without knowing the state space.

## Three algorithms available:
## social: https://github.com/WestleyArgentum/GeneticAlgorithms.jl
## with fitness based on parent's (so fewer MCs per evaluation are
## needed).
## biological: https://github.com/robertfeldt/BlackBoxOptim.jl
## sampled: Deterministic optimization under different Monte Carlo
## draws.

##### Social Genetics #####

using GeneticAlgorithms

module MimiGA

using GeneticAlgorithms
using Distributions

paramsize = 0
lowers = Array(Float64, 0)
uppers = Array(Float64, 0)
objective = (parameters) -> 0
iterations = 0
lifeseed = 0

function prepare(lows::Vector{Float64}, his::Vector{Float64}, obj::Function)
    global paramsize, lowers, uppers, objective, iterations

    paramsize = length(lows)
    lowers = lows
    uppers = his
    objective = obj
    iterations = 0
    lifeseed = 1
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
    srand(lifeseed)
    fit = objective(ent.parameters)
    if ent.fitness != nothing
        ent.fitness * .9 + fit * .1
    else
        fit * .1
    end
end

function group_entities(pop)
    global iterations, lifeseed

    iterations += 1
    if iterations % 100 == 0
        println(pop[1])
	if iterations >= 10000
	    return
	end
    end

    lifeseed += 1

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

##### Biological Genetics #####

using BlackBoxOptim

##########

type UncertainOptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    montecarlo::Function
    objective::Function
    lowers::Vector{Float64}
    uppers::Vector{Float64}
end

type UncertainOptimizationSolution
    problem::UncertainOptimizationProblem
    fmean::Float64
    fserr::Float64
    xmean::Vector{Float64}
    xserr::Vector{Float64}
    extra
end

"""Setup an optimization over Monte Carlo uncertainty."""
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)
    my_objective = unaryobjective(model, components, names, objective)

    UncertainOptimizationProblem(model, components, names, montecarlo, my_objective, my_lowers, my_uppers)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem, generator::Function, algorithm::Symbol, mcperlife::Int64, samples::Int64)
    sampleseed = 0

    function sample_objective(parameters::Vector{Float64})
        srand(sampleseed)

        total = 0
        for iter in 1:mcperlife
            optprob.montecarlo(optprob.model)
            total += optprob.objective(parameters)
        end

        total / mcperlife
    end

    if algorithm == :social
        MimiGA.prepare(optprob.lowers, optprob.uppers, sample_objective)
        gamodel = runga(MimiGA; initial_pop_size=samples)

        # get the latest population when the GA exited
        population = GeneticAlgorithms.population(gamodel)

        fmean = mean([entity.fitness for entity in population])
        fserr = std([entity.fitness for entity in population]) / sqrt(length(population))
        xmean = vec(mean(hcat([entity.parameters for entity in population]...), 2))
        xserr = vec(std(hcat([entity.parameters for entity in population]...), 2) / sqrt(length(population)))
        extra = population

    elseif algorithm == :biological
        res = bboptimize(p -> -sample_objective(p); SearchRange=collect(zip(optprob.lowers, optprob.uppers)), MaxFuncEvals=samples)

        fmean = NaN
        fserr = NaN
        xmean = best_candidate(res)
        xserr = repmat([NaN], length(optprob.lowers))
        extra = res
    elseif algorithm == :sampled
        opt = Opt(:LN_COBYLA, length(optprob.lowers))
        lower_bounds!(opt, optprob.lowers)
        upper_bounds!(opt, optprob.uppers)
        xtol_rel!(opt, minimum(1e-6 * (optprob.uppers - optprob.lowers)))
        max_objective!(opt, (p, g) -> sample_objective(p))

        sampleprob = OptimizationProblem(optprob.model, optprob.components, optprob.names, opt, Function[])

        allparams = Vector{Float64}[]

        for sample in 1:samples
            sampleseed = sample
            
            minf, minx = solution(sampleprob, generator, maxiter=10000)
            push!(allparams, Float64[minf; minx])
        end

        allparams = transpose(hcat(allparams...))

        fmean = mean(allparams[:, 1])
        fserr = std(allparams[:, 1]) / sqrt(samples)

        xmean = vec(mean(allparams[:, 2:end], 1))
        xserr = vec(std(allparams[:, 2:end], 1) / sqrt(samples))

        extra = nothing
    else
        error("Unknown algorithm")
    end

    UncertainOptimizationSolution(optprob, fmean, fserr, xmean, xserr, extra)
end

function evaluate(optprob::UncertainOptimizationProblem, param::Vector{Float64}, mcperlife::Int64, sampleseed=1)
    srand(sampleseed)

    total = 0
    for iter in 1:mcperlife
        optprob.montecarlo(optprob.model)
        total += optprob.objective(param)
    end

    total / mcperlife
end
