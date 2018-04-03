using Mimi
using OptiMimi
using Base.Test

@defcomp Simple begin
    region = Index()

    xx = Parameter(index=[region, time])
end

m = Model()
setindex(m, :region, collect(1:2))
setindex(m, :time, collect(1:3))

simple = addcomponent(m, Simple)
simple[:xx] = reshape(repeat([1, 2], outer=3), (2, 3))

gen(rr) = m.external_parameters[:xx].values[rr, 1]
hall = hallsingle(m, :Simple, :xx, gen, [:time])

@test reshape(hall.f, (2, 3)) == m.external_parameters[:xx].values
