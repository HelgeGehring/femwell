using Gridap
using Gridap.Geometry
using GridapGmsh
using GridapMakie, GLMakie
using PyCall

GLMakie.inline!(true)

py"
from collections import OrderedDict

import shapely

from femwell.mesh import mesh_from_OrderedDict

line1 = shapely.LineString(((0, 0), (0, 1)))
line2 = shapely.LineString(((7, 0), (7, 1)))
line_bottom = shapely.LineString(((-5, 0), (15, 0)))
line_top = shapely.LineString(((-5, 1), (15, 1)))
pml1 = shapely.box(-5, 0, 0, 1)
waveguide1 = shapely.box(0, 0, 5, 1)
waveguide2 = shapely.box(5, 0, 15, 1)
mesh_from_OrderedDict(
    OrderedDict(
        line1=line1,
        line2=line2,
        line_top=line_top,
        line_bottom=line_bottom,
        pml1=pml1,
        waveguide2=waveguide2,
        waveguide1=waveguide1,
    ),
    resolutions=dict(),
    default_resolution_max=0.05,
    filename='test.msh',
)
"

model = GmshDiscreteModel("test.msh")
Ω = Triangulation(model)
dΩ = Measure(Ω, 2)
Γ_in = BoundaryTriangulation(model, tags = "line1")
dΓ_in = Measure(Γ_in, 2)
Γ_out = BoundaryTriangulation(model, tags = "line2")
dΓ_out = Measure(Γ_out, 2)
#fig = plot(Ω)
#fig.axis.aspect = DataAspect()
#wireframe!(Ω, color = :black, linewidth = 1)
#display(fig)


mu_r = 1
k0 = 3.5

labels = get_face_labeling(model)
epsilons = ["waveguide1" => 1.0^2, "waveguide2" => 2^2, "pml1" => 1.0^2]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

h_m(y, b, m) = sqrt(((m == 0) ? 1 : 2) / b) * sin(m * π * y / b)
gamma_m(b, m) = sqrt(Complex((m * π / b)^2 - k0^2))

ε_total = (ε ∘ τ) + (x -> 0.1im * max(0, x[1] - 10)^2) + (x -> 0.1im * min(0, x[1] + 0.5)^2)

function lhs_maxwell(u, v)
    #println((h_m(get_cell_points(Ω)⋅VectorValue([0,1]), 1, 2)) ⋅ u)
    #println((x -> h_m(x[2], 1, 1)) ⋅ v)
    0.5 * ∫((1 / ε_total * ∇(u) ⋅ ∇(v)) - k0^2 * u ⋅ v)dΩ #+ 
    #∫(
    #    0.5 * sum(
    #        u * (x -> h_m(x[2], 1, m)) * gamma_m(1, m) * ((x -> h_m(x[2], 1, m)) ⋅ v) for
    #        m in range(1,1)
    #    ),
    #)dΓ_in
end

input_form_rhs(v) = ∫(-2 * gamma_m(1, 1) * ((x -> h_m(x[2], 1, 0) + 1) ⋅ v))dΓ_in

V = TestFESpace(
    model,
    ReferenceFE(lagrangian, Float64, 3),
    vector_type = Vector{ComplexF64},
    dirichlet_tags = Int[],# ["line_top", "line_bottom"],
)
U = TrialFESpace(V, x -> 0)

op = AffineFEOperator(lhs_maxwell, input_form_rhs, U, V)
sol = solve(op)

fig = plot(Ω, real(sol))
display(fig)
fig = plot(Ω, abs(sol * sol))
display(fig)

sum_in = [abs(sum(∫((x -> h_m(x[2], 1, m)) * sol)dΓ_in)) for m = 1:3]
sum_out = [abs(sum(∫((x -> h_m(x[2], 1, m)) * sol)dΓ_out)) for m = 1:3]
println(gamma_m(1, 1))
println(sum_in)
println(sum_out)

display(plot([abs(sol(Gridap.Point(x, 0.5))) for x = -5:0.02:15]))

println(
    (
        minimum([abs(sol(Gridap.Point(x, 0.5))) for x = 1:0.01:4]) + maximum([abs(sol(Gridap.Point(x, 0.5))) for x = 1:0.01:4])
    ) / 2,
)
println(minimum([abs(sol(Gridap.Point(x, 0.5))) for x = 1:0.01:4]))
println(maximum([abs(sol(Gridap.Point(x, 0.5))) for x = 1:0.01:4]))
println(maximum([abs(sol(Gridap.Point(x, 0.5))) for x = 6:0.01:8]))
