module Electrostatic

using Gridap
using Gridap.Algebra
using GridapPETSc

export compute_potential, current_density, power_density, potential

struct Potential
    potential::FEFunction
    conductivity::CellField
end

function compute_potential(
    conductivity::CellField,
    voltages::Dict{String,Float64};
    order::Int = 1,
    solver::LinearSolver = BackslashSolver(),
)
    tags = collect(keys(voltages))
    tags_voltages = [voltages[tag] for tag in tags]

    model = get_active_model(get_triangulation(conductivity))
    V0 = TestFESpace(
        model,
        ReferenceFE(lagrangian, Float64, order);
        conformity = :H1,
        dirichlet_tags = tags,
    )
    Ug = TrialFESpace(V0, tags_voltages)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, order)

    a(u, v) = ∫(conductivity * ∇(v) ⋅ ∇(u))dΩ
    b(v) = 0

    op = AffineFEOperator(a, b, Ug, V0)
    potential = solve(solver, op)

    return Potential(potential, conductivity)
end

potential(potential::Potential) = potential.potential
current_density(potential::Potential) =
    potential.conductivity * (norm ∘ ∇(potential.potential))
power_density(potential::Potential) =
    potential.conductivity *
    (norm ∘ ∇(potential.potential)) *
    (norm ∘ ∇(potential.potential))

end
