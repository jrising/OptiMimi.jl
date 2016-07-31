using Mimi
export grad_component

function grad_component(m::Model, component::Symbol, variable::Symbol, parameter::Symbol)
    fdf = unaryobjective(m, [component], [parameter], m -> vec(m[component, variable]))
    A = ForwardDiff.jacobian(fdf, convert(Vector{Float64}, vec(m.parameters[parameter].values)))
    LinearProgrammingRoom(component, variable, component, parameter, sparse(A))
end

function grad_component(indices::Dict{Symbol, Any}, addcomponent::Function, variable::Symbol, parameter::Symbol)
    m = Model(Number)
    for index in keys(indices)
        setindex(m, index, indices[index])
    end

    addcomponent(m, :Singleton)

    grad_component(m, :Singleton, variable, parameter)
end
