using Mimi
using OptiMimi
using Base.Test

# Create simple model
@defcomp Bellmano begin
    # The x-value to evaluate the quadratic
    consumption = Parameter(index=[time])

    # The y-value of the quadratic at the x-value
    utility = Variable(index=[index])
    disaster::Bool = Variable(index=[index])
end

function run_timestep(state::Bellmano, tt::Int64)
    v = state.Variables
    p = state.Parameters

    if tt == 1 || !v.disaster[tt-1]
        v.utility[tt] = sqrt(p.consumption[tt])
        v.disaster[tt] = rand() < p.consumption[tt]
    else
        v.utility[tt] = 0
        v.disaster[tt] = true
    end
end

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
