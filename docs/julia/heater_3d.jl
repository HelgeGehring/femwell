# ---
# jupyter:
#   jupytext:
#     custom_cell_magics: kql
#     formats: jl:percent,ipynb
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: base
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# # 3D thermal phase shifter

# %% tags=["remove-stderr"]
using CairoMakie

using Gridap
using GridapGmsh
using Gridap.Geometry
using GridapPETSc

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

dir = @__DIR__
run(`python $dir/heater_3d_mesh.py`)


model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)

labels = get_face_labeling(model)
tags = get_face_tag(labels, num_cell_dims(model))
τ = CellField(tags, Ω)


options = "-ksp_type cg -pc_type mg -ksp_monitor"
GridapPETSc.with(args = split(options)) do

    electrical_conductivity = [
        "core" => 1e-6,
        "box" => 1e-13,
        "clad" => 1e-13,
        "heater" => 2.3e6,
        "via2" => 3.77e7,
        "metal2" => 3.77e7,
        "via1" => 3.77e7,
        "metal3" => 3.77e7,
        "metal3#e1" => 3.77e7,
        "metal3#e2" => 3.77e7,
    ]
    electrical_conductivity =
        Dict(get_tag_from_name(labels, u) => v for (u, v) in electrical_conductivity)
    ϵ_electrical_conductivity(tag) = electrical_conductivity[tag]

    boundary_conditions = Dict(["metal3#e1___None" => 0.4, "metal3#e2___None" => 0.0])

    p0 = compute_potential(
        ϵ_electrical_conductivity ∘ τ,
        boundary_conditions,
        order = 1,
        solver = PETScLinearSolver(),
    )

    thermal_conductivities = [
        "core" => 90,
        "box" => 1.38,
        "clad" => 1.38,
        "heater" => 28,
        "via2" => 28,
        "metal2" => 28,
        "via1" => 28,
        "metal3" => 28,
        "metal3#e1" => 28,
        "metal3#e2" => 28,
    ]
    thermal_conductivities =
        Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_conductivities)
    ϵ_conductivities(tag) = thermal_conductivities[tag]

    temperatures = Dict("box___None" => 0.0)

    T0 = calculate_temperature(
        ϵ_conductivities ∘ τ,
        power_density(p0),
        temperatures,
        order = 1,
        solver = PETScLinearSolver(),
    )

    Ω_w = Triangulation(model, tags = "core")
    dΩ_w = Measure(Ω_w, 1)

    println("Average Temperature: ", ∑(∫(temperature(T0))dΩ_w) / ∑(∫(1)dΩ_w), " K")

    writevtk(
        Ω,
        "results",
        cellfields = [
            "potential" => potential(p0),
            "current" => current_density(p0),
            "temperature" => temperature(T0),
        ],
    )

    thermal_diffisitivities = [
        "core" => 90 / 711 / 2330,
        "box" => 1.38 / 709 / 2203,
        "clad" => 1.38 / 709 / 2203,
        "heater" => 28 / 598 / 5240,
        "via2" => 28 / 598 / 5240,
        "metal2" => 28 / 598 / 5240,
        "via1" => 28 / 598 / 5240,
        "metal3" => 28 / 598 / 5240,
        "metal3#e1" => 28 / 598 / 5240,
        "metal3#e2" => 28 / 598 / 5240,
    ]
    thermal_diffisitivities =
        Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_diffisitivities)
    ϵ_diffisitivities(tag) = thermal_diffisitivities[tag]

    uₕₜ = calculate_temperature_transient(
        ϵ_conductivities ∘ τ,
        ϵ_diffisitivities ∘ τ,
        power_density(p0),
        temperatures,
        temperature(T0) * 0,
        2e-6,
        2e-4,
        solver = PETScLinearSolver(),
    )

    stophere

    #createpvd("poisson_transient_solution") do pvd
    #    for (uₕ, t) in uₕₜ
    #        pvd[t] = createvtk(
    #            Ω,
    #            "poisson_transient_solution_$t" * ".vtu",
    #            cellfields = ["u" => uₕ],
    #        )
    #    end
    #end
    sums = [(t, ∑(∫(u)dΩ_w) / ∑(∫(1)dΩ_w)) for (u, t) in uₕₜ]

    figure = Figure()
    ax = Axis(figure[1, 1], ylabel = "Temperature / K", xlabel = "time / ms")

    t, s = getindex.(sums, 1), getindex.(sums, 2)
    lines!(ax, t * 1e3, s)
    display(figure)
end
