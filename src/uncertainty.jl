## We use Monte Carlo methods, because: 1. The size of the state space
## >> parameter space.  So just explore parameter space.  2. We cannot
## run from the last period reliably without knowing the state space.

## Two algorithms available:
## biological: https://github.com/robertfeldt/BlackBoxOptim.jl
## sampled: Deterministic optimization under different Monte Carlo
## draws.

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
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function=(model) -> nothing)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)
    my_objective = unaryobjective(model, components, names, objective)

    UncertainOptimizationProblem(model, components, names, montecarlo, my_objective, my_lowers, my_uppers)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem, generator::Function, algorithm::Symbol, mcperlife::Int64, samples::Int64)
    function sample_objective(parameters::Vector{Float64})
        total = 0
        for iter in 1:mcperlife
            optprob.montecarlo(optprob.model)
            total += optprob.objective(parameters)
        end

        total / mcperlife
    end

    if algorithm == :biological
        res = bboptimize(p -> -sample_objective(p); SearchRange=collect(zip(optprob.lowers, optprob.uppers)), MaxFuncEvals=samples, Method=:separable_nes)

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
