using Mimi
using OptiMimi
using Base.Test

# Basic model to optimize:
# Each period, you have 1 unit to divide between consuming and staving off disaster.
# If disaster happens, you get 0 consumption for the remainder of the game.
# Utility is defined as sqrt(c).  The probability of disaster is 1-c in each period.
# Assume 5% discounting and 10 periods.

## Analytical solution

## V_t = max_c u(c) + gamma (1 - c) V_t+1
## u'(c) - gamma V_t+1 = 0
## Let u(c) = sqrt(c), u'(c) = .5 c^(-.5)
## c = 1/(2 gamma V_t+1)^2
gamma = exp(-.05) # 5% discount rate

# consume everything in last period
VV_reverse = [1.]
consumption_reverse = [1.]

for tt in 10:-1:1
    consumption = 1 / (2 * gamma * VV_reverse[end])^2
    VV = sqrt(consumption) + gamma * (1 - consumption) * VV_reverse[end]
    push!(consumption_reverse, consumption)
    push!(VV_reverse, VV)
end

println(reverse(consumption_reverse))
println(reverse(VV_reverse))

# OptiMimi solution

# Create simple model
@defcomp Bellmano begin
    # The x-value to evaluate the quadratic
    consumption = Parameter(index=[time])

    # The y-value of the quadratic at the x-value
    utility = Variable(index=[time])
    disaster::Bool = Variable(index=[time])
    bonus = Variable()
end

function run_timestep(state::Bellmano, tt::Int64)
    v = state.Variables
    p = state.Parameters

    v.bonus = 0
    if tt == 1 || !v.disaster[tt-1]
        v.utility[tt] = sqrt(p.consumption[tt])
        v.disaster[tt] = rand() < p.consumption[tt]
	if tt == 10 && !v.disaster[tt]
	    v.bonus = 1
	end
    else
        v.utility[tt] = 0
        v.disaster[tt] = true
    end
end

m = Model()
setindex(m, :time, collect(1:10))

bellmano = addcomponent(m, Bellmano)
bellmano[:consumption] = repmat([.25], 10)

run(m)
m[:Bellmano, :utility]

import OptiMimi.uncertainproblem
using DataFrames

df = DataFrame(mcperlife=[], cons1=[], cons2=[], cons3=[], cons4=[], cons5=[], cons6=[], cons7=[], cons8=[], cons9=[], cons10=[])
push!(df, [Inf; reverse(consumption_reverse)[1:10]])

for mcperlife in 1:5:40
    println(mcperlife)
    prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(m[:Bellmano, :utility] .* exp(-(0:9) * .05)) + m[:Bellmano, :bonus] * exp(-10 * .05), (model) -> nothing, mcperlife)
    soln = solution(prob)
    push!(df, [mcperlife; soln])
    println(df)
end

writetable("attempts.csv", df)

for mcperlife in 50:10:100
    println(mcperlife)
    prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(m[:Bellmano, :utility] .* exp(-(0:9) * .05)) + m[:Bellmano, :bonus] * exp(-10 * .05), (model) -> nothing, mcperlife)
    soln = solution(prob)
    push!(df, [mcperlife; soln])
    println(df)
end

writetable("attempts.csv", df)

for mcperlife in 120:20:200
    println(mcperlife)
    prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(m[:Bellmano, :utility] .* exp(-(0:9) * .05)) + m[:Bellmano, :bonus] * exp(-10 * .05), (model) -> nothing, mcperlife)
    soln = solution(prob)
    push!(df, [mcperlife; soln])
    println(df)
end

writetable("attempts.csv", df)

using MultivariateStats

llsq(convert(Matrix{Float64}, [ones(nrow(df)-1) log(convert(Vector{Float64}, df[:mcperlife][2:end]))]), convert(Vector{Float64}, df[:cons1][2:end]); bias=false)

finite = df[:cons1] .> df[1, :cons1]

llsq(convert(Matrix{Float64}, [ones(sum(finite)) convert(Vector{Float64}, df[:mcperlife][finite])]), convert(Vector{Float64}, log(df[:cons1][finite] - df[1, :cons1])); bias=false)
