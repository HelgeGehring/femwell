using Gridap
using Gridap.Geometry
using GridapGmsh
using GridapMakie, GLMakie
using PyCall

GLMakie.inline!(true)

solver = LinearFESolver(BackslashSolver())
solver = LinearFESolver(LUSolver())

model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)
dΩ = Measure(Ω, 4)
Γ_in = BoundaryTriangulation(model, tags = "cladding___PML_cladding")
dΓ_in = Measure(Γ_in, 4)
Γ_out = BoundaryTriangulation(model, tags = "substrate___PML_substrate")
dΓ_out = Measure(Γ_out, 4)
#fig = plot(Ω)
#fig.axis.aspect = DataAspect()
#wireframe!(Ω, color = :black, linewidth = 1)
#display(fig)


mu_r = 1
k0 = 10

labels = get_face_labeling(model)

n1 = 1
n2 = 1.2
epsilons = [
    "substrate" => n1^2,
    "cladding" => n2^2,
    "PML_substrate" => n1^2,
    "PML_cladding" => n2^2,
]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

ε_total = (ε ∘ τ) + (x -> 10im * max(0, x[3] - 1)^2) + (x -> 10im * min(0, x[3] + 1)^2)

function lhs_maxwell(u, v)
    0.5 * ∫((1 / ε_total * curl(u) ⋅ curl(v)) - k0^2 * u ⋅ v)dΩ
end

input_form_rhs(v) = ∫(-2 * ((x -> VectorValue(1, 0, 0)) ⋅ v))dΓ_in

V = TestFESpace(
    model,
    ReferenceFE(nedelec, Float64, 1),
    conformity = :HCurl,
    vector_type = Vector{ComplexF64},
    #dirichlet_tags = ["cladding___boundary", "substrate___boundary"],
    #dirichlet_masks=[(true,false,false),(true,false,false)]
)
U = TrialFESpace(V, x -> VectorValue(0, 0, 0))
println(num_free_dofs(V))

op = AffineFEOperator(lhs_maxwell, input_form_rhs, U, V)

sol = solve(solver, op)

Base.real(x::VectorValue{3,ComplexF64}) = VectorValue(real.(x.data))
Base.imag(x::VectorValue{3,ComplexF64}) = VectorValue(imag.(x.data))
writevtk(
    Ω,
    "results",
    cellfields = [
        "epsilon_real" => ε ∘ τ,
        "epsilon_imag" => imag(ε_total),
        "field_real" => real(sol),
        "field_imag" => imag(sol),
    ],
)

#fig = plot(Ω, real(sol))
#display(fig)
#fig = plot(Ω, abs(sol ⋅ sol))
#display(fig)

sum_in = abs(sum(∫(norm ∘ sol)dΓ_in))
sum_out = abs(sum(∫(norm ∘ sol)dΓ_out))
println(sum_in)
println(sum_out)

display(
    plot([imag(sol(Gridap.Point(0.5, 0.5, x)) ⋅ VectorValue(1, 0, 0)) for x = -2:0.01:2]),
)

println(
    (
        minimum([abs((norm ∘ sol)(Gridap.Point(0.5, 0.5, x))) for x = -1:0.01:1]) + maximum([abs((norm ∘ sol)(Gridap.Point(0.5, 0.5, x))) for x = -1:0.01:1])
    ) / 2,
)
println(minimum([abs((norm ∘ sol)(Gridap.Point(0.5, 0.5, x))) for x = -1:0.01:1]))
println(maximum([abs((norm ∘ sol)(Gridap.Point(0.5, 0.5, x))) for x = -1:0.01:1]))
println(maximum([abs((norm ∘ sol)(Gridap.Point(0.5, 0.5, x))) for x = -1:0.01:1]))

transmission =
    maximum([imag(sol(Gridap.Point(0.5, 0.5, x)) ⋅ VectorValue(1, 0, 0)) for x = -1:0.01:0])
in_plus_reflection =
    maximum([imag(sol(Gridap.Point(0.5, 0.5, x)) ⋅ VectorValue(1, 0, 0)) for x = 0:0.01:1])

println(transmission / in_plus_reflection)

reflection = abs((n1 - n2) / (n1 + n2))
println((1 - reflection) / (1 + reflection))
