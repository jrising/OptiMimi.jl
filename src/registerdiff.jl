type RegisterGradient
    libeval::Function

    dx::Float64
    funcs::Vector{Function}
end

function unarylibrary(model::Model, components::Vector{Symbol}, names::Vector{Symbol})
    function evaluator(xx::Vector, funcs::Vector{Function})
        setparameters(model, components, names, xx)
        run(model)

        map(func -> func(model), funcs)
    end

    evaluator
end


RegisterGradient(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, dx::Float64) =
    RegisterGradient(unarylibrary(model, components, names), dx, Function[])

function addfunction!(rg::RegisterGradient, func::Function)
    rg.funcs = [rg.funcs; func]
end

function getgradients(rg::RegisterGradient, xx::Vector{Float64})
    baseline = rg.libeval(xx, rg.funcs)

    gradient = Array{Float64}(length(rg.funcs), length(xx))

    for ii in 1:length(xx)
        println(ii / length(xx))
        # Adjust a single value
        xx[ii] = xx[ii] + rg.dx
        gradient[:, ii] = (rg.libeval(xx, rg.funcs) - baseline) / rg.dx
        xx[ii] = xx[ii] - rg.dx
    end

    baseline, gradient
end
