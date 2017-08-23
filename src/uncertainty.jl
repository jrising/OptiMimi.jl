type UncertainOptimizationProblem
    model::Model
    components::Vector{Symbol}
    names::Vector{Symbol}
    montecarlo::Function
    objective::Function
    lowers::Vector{Float64}
    uppers::Vector{Float64}
end

"""Setup an optimization over Monte Carlo uncertainty."""
function uncertainproblem(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, montecarlo::Function, mcperlife::Int64)
    my_lowers, my_uppers, totalvars = expandlimits(model, components, names, lowers, uppers)
    my_objective = unaryobjective(model, components, names, objective)

    UncertainOptimizationProblem(model, components, names, montecarlo, my_objective, my_lowers, my_uppers)
end

"""Solve an optimization over Monte Carlo uncertainty."""
function solution(optprob::UncertainOptimizationProblem, generator::Function, samples::Int64)
    sampleseed = 0
    
    function sample_objective(parameters::Vector{Float64})
        srand(sampleseed)
        
        total = 0
        for iter in 1:mcperlife
            optprob.montecarlo(model)
            total += optprob.objective(parameters)
        end

        -total / mcperlife
    end

    opt = Opt(length(optprob.lowers), :LN_COBYLA)
    lower_bounds!(opt, optprob.lowers)
    upper_bounds!(opt, optprob.uppers)
    xtol_rel!(opt, minimum(1e-6 * (optprob.uppers - optprob.lowers)))
    max_objective!(opt, sample_objective)

    sampleprob = OptimizationProblem(optprob.model, optprob.components, optprob.names, opt, Function[])
    
    allparams = []

    for sample in 1:samples
        sampleseed = sample
        
        minf, minx = solution(sampleprob, generator, maxiter=10000)
        push(allparams, [minf; minx])
        println(allparams)
    end
end

