using Mimi
using OptiMimi
using Base.Test

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

    if rand() < p.consumption[tt]
        v.utility[tt] = 0
    else
        v.utility[tt] = sqrt(p.consumption[tt])
    end
end

m = Model()
setindex(m, :time, collect(1:10))

bellmano = addcomponent(m, Bellmano)
bellmano[:consumption] = repmat([.33], 10) # XXX: Try initializing at the right answer

run(m)
m[:Bellmano, :utility]

import OptiMimi.uncertainproblem
using DataFrames

df = DataFrame(mcperlife=[], cons1=[], cons2=[], cons3=[], cons4=[], cons5=[], cons6=[], cons7=[], cons8=[], cons9=[], cons10=[])
push!(df, [Inf; repmat([1/3], 10)])

for mcperlife in 1:5:41
    println(mcperlife)
    prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(sqrt(m[:Bellmano, :utility]) .* exp(-(0:9) * .05)), (model) -> nothing, mcperlife)
    soln = solution(prob)
    push!(df, [mcperlife; soln])
    println(df)
end

writetable("attempts.csv", df)

for mcperlife in 50:10:100
    println(mcperlife)
    prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(sqrt(m[:Bellmano, :utility]) .* exp(-(0:9) * .05)), (model) -> nothing, mcperlife)
    soln = solution(prob)
    push!(df, [mcperlife; soln])
    println(df)
end

writetable("attempts.csv", df)

for mcperlife in 120:20:200
    println(mcperlife)
    prob = uncertainproblem(m, [:Bellmano], [:consumption], [0.], [1.], m -> sum(sqrt(m[:Bellmano, :utility]) .* exp(-(0:9) * .05)), (model) -> nothing, mcperlife)
    soln = solution(prob)
    push!(df, [mcperlife; soln])
    println(df)
end

writetable("attempts.csv", df)

using MultivariateStats

llsq(convert(Matrix{Float64}, [ones(nrow(df)-1) log(convert(Vector{Float64}, df[:mcperlife][2:end]))]), convert(Vector{Float64}, df[:cons1][2:end]); bias=false)

finite = df[:cons1] .> df[1, :cons1]

llsq(convert(Matrix{Float64}, [ones(sum(finite)) convert(Vector{Float64}, df[:mcperlife][finite])]), convert(Vector{Float64}, log(convert(Vector{Float64}, df[:cons1][finite]) - df[1, :cons1])); bias=false)
