module Electrostatic

using Gridap
using Gridap.Algebra
using GridapMakie, CairoMakie

export compute_potential, current_density, power_density, potential, plot_potential

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

function plot_potential(potential::Potential)
    fig, _, plt = plot(get_triangulation(potential.potential), potential.potential)
    Colorbar(fig[1, 2], plt)
    display(fig)
end

end
