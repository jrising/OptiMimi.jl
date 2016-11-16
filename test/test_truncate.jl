using Mimi
using OptiMimi
using Base.Test

import OptiMimi.gettruncdims, OptiMimi.gettruncdimnames, OptiMimi.untruncindex

@defcomp Example begin
    one = Index()
    two = Index()

    data = Parameter(index=[time, one, two])
    output = Parameter(index=[time])
end

m = Model()
setindex(m, :time, collect(1:2))
setindex(m, :one, collect(1:3))
setindex(m, :two, collect(1:4))

example = addcomponent(m, Example)
example[:data] = ones(2, 3, 4)

@test gettruncdims(m, :Example, :data, [:one]) == [2, 4]
@test gettruncdims(m, :Example, :data, [:time]) == [3, 4]

@test gettruncdimnames(m, :Example, :data, [:one]) == [:time, :two]
@test gettruncdimnames(m, :Example, :data, [:time]) == [:one, :two]

@test untruncindex(1, m, :Example, :data, [:two]) == (1, 1, :)
@test untruncindex(2, m, :Example, :data, [:two]) == (2, 1, :)
@test untruncindex(3, m, :Example, :data, [:two]) == (1, 2, :)

hall = hallsingle(m, :Example, :data, (tt, ii, jj) -> tt * ii)
@test trunchall_dropdim(m, hall, :two, sum).f == [1*1 * 4, 2*1 * 4, 1*2 * 4, 2*2 * 4, 1*3 * 4, 2*3 * 4]
@test trunchall_dropdim(m, hall, :two, reduce_same).f == [1*1, 2*1, 1*2, 2*2, 1*3, 2*3]
