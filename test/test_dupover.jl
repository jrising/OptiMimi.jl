using Mimi
using OptiMimi
using Test

@defcomp Simple begin
    region = Index()

    xx = Parameter(index=[time, region])
end

m = Model()
set_dimension!(m, :region, collect(1:2))
set_dimension!(m, :time, collect(1:3))

simple = add_comp!(m, Simple)
xx = reshape(repeat([1, 2], inner=3), (3, 2))
simple[:xx] = xx

gen(rr) = m.md.external_params[:xx].values[rr, 1]
hall = hallsingle(m, :Simple, :xx, gen, [:time])

@test reshape(hall.f, (3, 2)) == m.md.external_params[:xx].values
