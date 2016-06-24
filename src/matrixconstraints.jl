export savelpconstraints, MatrixConstraintSet

type MatrixConstraintSet
    components::Vector{Symbol}
    names::Vector{Symbol}
    lowers::Vector{Float64}
    uppers::Vector{Float64}
    f::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int64}
    b::Vector{Float64}
end

function lpconstraints(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, objective::Function, constraints::Vector{Function}; verbose=false)
    if verbose
        println("Making zero point.")
    end
    initial = make0(model, names)

    if model.numberType == Number
        myobjective = unaryobjective(model, components, names, objective)
        f = ForwardDiff.gradient(myobjective, initial)

        myconstraints = Function[]
        for constraint in constraints
            myconstraints = [myconstraints; unaryobjective(model, components, names, constraint)]
        end

        # Copied from makematrix.jl
        # Could be made more efficient to get b values from gradient
        A = Float64[]
        b = Float64[]
        for myconstraint in myconstraints
            A = [A; ForwardDiff.gradient(myconstraint, initial)]
            b = [b; myconstraint(initial)] # Not the final values for 'b' yet!
        end

        A = reshape(A, (length(initial), div(length(A), length(initial))))'
    else
        if verbose
            println("Adding to gradient collection.")
        end
        rg = RegisterGradient(model, components, names, 1.0)
        addfunction!(rg, objective)
        for constraint in constraints
            addfunction!(rg, constraint)
        end

        vb, fA = getgradients(rg, initial, verbose=verbose, usesparse=false)
        f = vec(fA[1, :])
        A = fA[2:end, :]
        b = vb[2:end]
    end

    b = A * initial - b

    f, b, A
end

function savelpconstraints(model::Model, components::Vector{Symbol}, names::Vector{Symbol}, lowers::Vector{Float64}, uppers::Vector{Float64}, objective::Function, constraints::Vector{Function}=Function[]; verbose=false)
    f, b, A = lpconstraints(model, components, names, objective, constraints, verbose=verbose)

    sA = sparse(A)
    MatrixConstraintSet(components, names, lowers, uppers, f, sA, b)
end

function combineconstraints(f::Vector{Float64}, b::Vector{Float64}, A::Matrix{Float64}, model::Model, components::Vector{Symbol}, names::Vector{Symbol}, constraints::Vector{MatrixConstraintSet})
    for constraint in constraints
        f = [f; constraint.f] # these need to be variables never seen before!
        b = [b; constraint.b]
        A, components, names = addnamedcols(A, model, components, names, constraint.components, constraint.names)
        A = addnamedrows(A, model, components, names, constraint.components, constraint.names, full(constraint.A))
    end

    f, b, A, components, names
end

function addnamedcols(A::Matrix{Float64}, model::Model, components::Vector{Symbol}, names::Vector{Symbol}, pluscomponents::Vector{Symbol}, plusnames::Vector{Symbol})
    found = Bool[]
    for ii in 1:length(pluscomponents)
        found = [found; any((components == pluscomponents[ii]) & (names == plusnames[ii]))]
    end

    # Add 0's for any unfound
    for (ii, len, isscalar) in @task nameindexes(model, plusnames[!found])
        A = [A zeros(size(A)[1], len)]
    end

    A, [components; pluscomponents[!found]], [names; plusnames[!found]]
end

function addnamedrows(A::Matrix{Float64}, model::Model, components::Vector{Symbol}, names::Vector{Symbol}, pluscomponents::Vector{Symbol}, plusnames::Vector{Symbol}, plusA::Matrix{Float64})
    rows = zeros(size(plusA)[1], size(A)[2])

    # Get the indexes into plusA
    plusAindexes = Int64[]
    startindex = 1
    for (ii, len, isscalar) in @task nameindexes(model, plusnames)
        plusAindexes = [plusAindexes; startindex]
        startindex += len
    end
    plusAindexes = [plusAindexes; startindex]

    # Filter these into the right spots in rows
    startindex = 1
    for (ii, len, isscalar) in @task nameindexes(model, names)
        if names[ii] in plusnames
            which = find(plusnames .== names[ii])[1]
            println(startindex:startindex+len-1)
            println(plusAindexes[which])
            println(plusAindexes[which+1])
            rows[:, startindex:startindex+len-1] = plusA[:, plusAindexes[which]:plusAindexes[which+1]-1]
        end
    end

    [A; rows]
end

function combinelimits(exlowers::Vector{Float64}, exuppers::Vector{Float64}, model::Model, components::Vector{Symbol}, names::Vector{Symbol}, constraints::Vector{MatrixConstraintSet})
    for constraint in constraints
        for (ii, len, isscalar) in @task nameindexes(model, constraint.names)
            append!(exlowers, [constraint.lowers[ii] for jj in 1:len])
            append!(exuppers, [constraint.uppers[ii] for jj in 1:len])
        end
    end

    exlowers, exuppers
end
