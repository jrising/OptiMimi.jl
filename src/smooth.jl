using Interpolations
import Dierckx
using Distributions

using NLopt
# import OptiMimi.OptimizationProblem, OptiMimi.getdims, OptiMimi.gradfreeobjective
# using Plots

function piecewisesolution_single(optprob::OptimizationProblem, initgen::Function, maxiter=Inf, verbose=false; reltol=1e-6)
    len = getdims(optprob.model, optprob.components[1], optprob.names[1])[1]
    numsegs = Int64.([2 .^(0:(floor(log2(len-1))-1)); len-1])

    last_xx = last_soln = nothing

    for numseg in numsegs
        print(numseg)
        xx = trunc.(Int, range(1, stop=len, length=numseg+1))

        keeps = [true for ii in 1:length(xx)]
        if numseg > 2
            # Drop xx where [LinearInterp - CubicInterp] < reltol
            linears = LinearInterpolation(last_xx, last_soln)(xx)
            cubics = Dierckx.evaluate(Dierckx.Spline1D(last_xx, last_soln, k=min(3, length(last_xx) - 1)), xx)
            keeps = abs.(linears .- cubics) ./ ((linears .+ cubics) / 2) .> reltol
            keeps[1] = true
            keeps[length(keeps)] = true

            # Drop knots only if can be linearly guessed
            for ii in 2:numseg
                if xx[ii] ∈ last_xx
                    jj = findfirst(last_xx .== xx[ii])
                    guess = LinearInterpolation([last_xx[jj-1], last_xx[jj+1]], [last_soln[jj-1], last_soln[jj+1]])(xx[ii])
                    if abs(guess - last_soln[jj]) / last_soln[jj] > reltol
                        keeps[ii] = true
                    end
                end
            end

            println(keeps)
        end

        paramtrans = convert(Vector{Function}, [yy -> LinearInterpolation(xx[keeps], yy)(1:len)])

        opt = Opt(optprob.opt.algorithm, sum(keeps))
        lower_bounds!(opt, optprob.opt.lower_bounds[xx[keeps]])
        upper_bounds!(opt, optprob.opt.upper_bounds[xx[keeps]])
        xtol_rel!(opt, reltol)

        myobjective = gradfreeobjective(optprob.model, optprob.components, optprob.names, optprob.objective, paramtrans)

        max_objective!(opt, myobjective)

        for constraint in optprob.constraints
            let this_constraint = constraint
                function my_constraint(xx::Vector, grad::Vector)
                    setparameters(optprob.model, optprob.components, optprob.names, xx, paramtrans)
                    run(optprob.model)
                    this_constraint(optprob.model)
                end

                inequality_constraint!(opt, my_constraint)
            end
        end

        if sum(keeps) == len - 1
            subprob = OptimizationProblem(optprob.model, optprob.components, optprob.names, optprob.objective, opt, optprob.constraints, nothing)
        else
            subprob = OptimizationProblem(optprob.model, optprob.components, optprob.names, optprob.objective, opt, optprob.constraints, paramtrans)
        end

        if isnothing(last_xx)
            lin_yy = nothing
        else
            lin_yy = LinearInterpolation(last_xx, last_soln)(xx[keeps])
        end
        soln = solution(subprob, initgen(sum(keeps), lin_yy))

        if numseg > 2
            last_soln = LinearInterpolation(last_xx, last_soln)(xx)
            last_soln[keeps] = soln[2]
        else
            last_soln = soln[2]
        end
        last_xx = xx
        println([last_xx, last_soln])
    end

    last_soln
end
