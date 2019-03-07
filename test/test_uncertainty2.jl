using Mimi
using OptiMimi
using Test

# Easier basic model to optimize:
# Each period, you have 1 unit to divide between consuming and staving off disaster.
# If disaster happens, you get 0 consumption for that period.
# Utility is defined as sqrt(c).  The probability of disaster is 1-c in each period.
# Assume 5% discounting and 10 periods.

## Analytical solution in 1 period

## max_c (1 - c) u(c) => c = 1/3
gamma = exp(-.05) # 5% discount rate

# OptiMimi solution

# Create simple model
@defcomp Bellmano begin
    # The x-value to evaluate the quadratic
    consumption = Parameter(index=[time])

    # The y-value of the quadratic at the x-value
    utility = Variable(index=[time])
end

function run_timestep(state::Bellmano, tt::Int64)
    v = state.Variables
    p = state.Parameters

    #v.utility[tt] = sqrt(p.consumption[tt]) * (1 - p.consumption[tt])

    if rand() < p.consumption[tt]
        v.utility[tt] = 0
    else
        v.utility[tt] = sqrt(max(0, p.consumption[tt]))
    end
end

m = Model()
setindex(m, :time, collect(1:10))

bellmano = addcomponent(m, Bellmano)
bellmano[:consumption] = repmat([.25], 10)

run(m)

import OptiMimi.uncertainproblem
using DataFrames

df = DataFrame(algorithm=[], mcperlife=[], cons1=[], cons2=[], cons3=[], cons4=[], cons5=[], cons6=[], cons7=[], cons8=[], cons9=[], cons10=[], cons1serr=[], cons2serr=[], cons3serr=[], cons4serr=[], cons5serr=[], cons6serr=[], cons7serr=[], cons8serr=[], cons9serr=[], cons10serr=[], time=[])
push!(df, [:exact; Inf; repmat([1/3], 10); repmat([0], 10); 0])

for samplemult in [1, 2, 4]
    for mcperlife in [5, 20, 50]
        println(mcperlife)
        prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(m[:Bellmano, :utility] .* exp(-(0:9) * .05)), (model) -> nothing)

        tic()
        soln = solution(prob, () -> repmat([.25], 10), :social, mcperlife, samplemult)
        time = toc()
        push!(df, [:social; mcperlife; soln.xmean; soln.xserr; time])

        println(df)

        tic()
        soln = solution(prob, () -> repmat([.25], 10), :biological, mcperlife, 120000*samplemult)
        time = toc()
        push!(df, [:biological; mcperlife; soln.xmean; soln.xserr; time])

        println(df)

        tic()
        soln = solution(prob, () -> repmat([.25], 10), :sampled, mcperlife, 300*samplemult)
        time = toc()
        push!(df, [:sampled; mcperlife; soln.xmean; soln.xserr; time])

        println(df)
    end

    writetable("attempts-compare2.csv", df)
end
