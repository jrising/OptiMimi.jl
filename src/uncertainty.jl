## We use Monte Carlo methods, because: 1. The size of the state space
## >> parameter space.  So just explore parameter space.  2. We cannot
## run from the last period reliably without knowing the state space.

## Two algorithms available:
## biological: https://github.com/robertfeldt/BlackBoxOptim.jl
## sampled: Deterministic optimization under different Monte Carlos
## montecarlo: Optimization of Monte Carlo mean of objective
## draws.

##### Biological Genetics #####

import BlackBoxOptim
using Random

##########

mutable struct UncertainOptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    montecarlo::Union{Function, SimulationDef}
    objective::Function
    lowers::Vector{Float64}
    uppers::Vector{Float64}
end

mutable struct UncertainOptimizationSolution
    problem::UncertainOptimizationProblem
    fmean::Float64
    fserr::Float64
    xmean::Vector{Float64}
    xserr::Vector{Float64}
    extra
end

"""Setup an optimization over Monte Carlo uncertainty."""
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function=(model) -> nothing)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)
    UncertainOptimizationProblem(model, components, names, montecarlo, objective, my_lowers, my_uppers)
end

function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, sim::SimulationDef, mcnum=Int64)
    function objectivetopayload(sim_inst::SimulationInstance, trialnum::Int, ntimesteps::Int, tup::Nothing)
        model = sim_inst.models[1]
        results = Mimi.payload(sim_inst)
        results[trialnum] = objective(model)
    end

    function uncertainobjective(model::Model)
        Mimi.set_payload!(sim, zeros(mcnum))
        Random.seed!(0)
        mcout = run(sim, model, mcnum, post_trial_func=objectivetopayload)
        values = Mimi.payload(mcout)
        mean(values[isfinite.(values)])
    end

    problem(model, components, names, lowers, uppers, uncertainobjective)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem, generator::Function, algorithm::Symbol, mcperlife::Int64, samples::Int64)
    if algorithm ∈ [:biological, :sampled]
        my_objective = unaryobjective(model, components, names, objective)

        function sample_objective(parameters::Vector{Float64})
            total = 0
            for iter in 1:mcperlife
                optprob.montecarlo(optprob.model)
                total += optprob.objective(parameters)
            end

            total / mcperlife
        end

        if algorithm == :biological
            res = BlackBoxOptim.bboptimize(p -> -sample_objective(p); SearchRange=collect(zip(optprob.lowers, optprob.uppers)), MaxFuncEvals=samples, Method=:separable_nes)

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

            sampleseed = 0
            function hold_objective(parameters::Vector{Float64}, grad::Vector{Float64})
                srand(sampleseed)
                sample_objective(parameters)
            end

            max_objective!(opt, hold_objective)

            sampleprob = OptimizationProblem(optprob.model, optprob.components, optprob.names, optprob.objective, opt, Function[])

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
        end
    elseif algorithm == :montecarlo
        minf, minx = solution(subprob, generator, samples)
        fmean = minf
        fserr = 0
        xmean = minx
        xserr = 0
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
