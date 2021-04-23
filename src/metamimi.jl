"Return a vector of the indices defining the parameter or variable."
function getdims(model::Model, component::Symbol, name::Symbol)
    dims = getdimnames(model, component, name)
    Vector{Int64}(map(dim -> dim_count(model, dim), dims))
end

"""
Return the symbols representing each of the dimensions for this variable or parameter.
"""
function getdimnames(model::Model, component::Symbol, name::Symbol)
    if name in parameter_names(model, component)
        parameter_dimensions(model, component, name)
    else
        variable_dimensions(model, component, name)
    end
end

"""
Avoid the distinction between set_param! and update_param!
"""
function set_or_update_param!(model::Model, component::Symbol, name::Symbol, value)
    try
        update_param!(model, name, value)
    catch
        set_param!(model, component, name, value)
    end
end

