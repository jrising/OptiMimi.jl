using ForwardDiff
using Mimi
export grad_component

# Optimize a univariate objective, by constructing a close linear approximation
function univariatematrix(func, lb, ub, x0, tolerance)
    # Collect all derivatives
    result = DerivativeResult(x0)
    ForwardDiff.derivative!(result, func, x0);

    pointslope = ForwardDiff.derivative(result)
    pointvalue = ForwardDiff.value(result)

    pointsecond = ForwardDiff.derivative(x -> ForwardDiff.derivative(func, x), x0)

    # Determine span for which linear is good approx
    pointsafespan = sqrt(2 * tolerance / abs(pointsecond))
    leftsafespan = min(pointsafespan, tolerance / abs(func(lb) - -pointslope * (x0 - lb)))
    rightsafespan = min(pointsafespan, tolerance / abs(func(ub) - pointslope * (ub - x0)))

    if lb < x0 - leftsafespan
        # Create average slope from here to edge
        leftslope = (func(x0) - func(x0 - leftsafespan)) / leftsafespan
        # Recurse for areas left of here
        leftfs, leftmins, leftmaxs = univariatematrix(func, lb, x0 - leftsafespan,
                                                      x0 - leftsafespan, tolerance)
    else
        leftslope = (func(x0) - func(lb)) / (x0 - lb)
        leftfs = leftmins = leftmaxs = []
    end

    if ub > x0 + rightsafespan
        # Create average slope from here to edge
        rightslope = (func(x0 + rightsafespan) - func(x0)) / rightsafespan
        # Recurse for areas right of here
        rightfs, rightmins, rightmaxs = univariatematrix(func, x0 + rightsafespan, ub,
                                                         x0 + rightsafespan, tolerance)
    else
        rightslope = (func(ub) - func(x0)) / (ub - x0)
        rightfs = rightmins = rightmaxs = []
    end

    if x0 == lb
        fs = [rightslope; leftfs; rightfs]
        mins = [0.; leftmins; rightmins]
        maxs = [min(rightsafespan, ub - x0); leftmaxs; rightmaxs]
    elseif x0 == ub
        fs = [-leftslope; leftfs; rightfs]
        mins = [0.; leftmins; rightmins]
        maxs = [min(leftsafespan, x0 - lb); leftmaxs; rightmaxs]
    else
        # leftslope is negative, because if upward sloping then this would take us down
        fs = [-leftslope; rightslope; leftfs; rightfs]
        mins = [0.; 0.; leftmins; rightmins]
        maxs = [min(leftsafespan, x0 - lb); min(rightsafespan, ub - x0); leftmaxs; rightmaxs]
    end

    fs, mins, maxs
end

function grad_component(m::Model, component::Symbol, variable::Symbol, parameter::Symbol)
    fdf = unaryobjective(m, [component], [parameter], m -> vec(m[component, variable]))
    A = ForwardDiff.jacobian(fdf, convert(Vector{Float64}, vec(m.external_parameters[parameter].values)))
    LinearProgrammingRoom(component, variable, component, parameter, sparse(A))
end

function grad_component(indices::Dict{Symbol, Any}, add_comp!::Function, variable::Symbol, parameter::Symbol)
    m = Model(Real)
    for index in keys(indices)
        set_dimension!(m, index, indices[index])
    end

    add_comp!(m, :Singleton)
    run(m) # Needs to be run for ModelInstance

    grad_component(m, :Singleton, variable, parameter)
end
