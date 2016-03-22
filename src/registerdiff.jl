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

function getgradients(rg::RegisterGradient, xx::Vector{Float64}; verbose=false, usesparse=true)
    baseline = rg.libeval(xx, rg.funcs)
    if any(isnan(baseline))
        error("Cannot compute the baseline for gradients.")
    end

    if usesparse
        gradients = spzeros(length(rg.funcs), length(xx))
        println(issparse(gradients))
        println(nnz(gradients))
    else
        gradients = Array{Float64}(length(rg.funcs), length(xx))
    end

    for ii in 1:length(xx)
        println(ii / length(xx))
        # Adjust a single value
        xx[ii] = xx[ii] + rg.dx

        newvals = rg.libeval(xx, rg.funcs)
        if any(isnan(newvals))
            error("Cannot calculate the gradient for constraint $ii")
        end

        if usesparse
            for jj in 1:length(rg.funcs)
                if newvals[jj] != baseline[jj]
                    gradients[jj, ii] = (newvals[jj] - baseline[jj]) / rg.dx
                end
            end
        else
            gradients[:, ii] = (newvals - baseline) / rg.dx
        end

        xx[ii] = xx[ii] - rg.dx
    end

    baseline, gradients
end
