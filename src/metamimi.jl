"Return a vector of the indices defining the parameter or variable."
function getdims(model::Model, component::Symbol, name::Symbol)
    indexes = getindexlabels(model, component, name)
    Vector{Int64}(map(index -> getindexcount(model, index), indexes))
end
