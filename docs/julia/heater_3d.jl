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

# %% tags=["remove-stderr", "hide-input", "hide-output"]
using CairoMakie
using Printf

using Gridap
using GridapGmsh
using Gridap.Geometry
using GridapPETSc

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

dir = @__DIR__
run(`python $dir/heater_3d_mesh.py`)

# %% [markdown]
# Let's start with defining the constants for the simulation:

# %% tags=["hide-output"]
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
thermal_conductivities = [
    "core" => 90,
    "box" => 1.38,
    "clad" => 1.38,
    "heater" => 400,
    "via2" => 400,
    "metal2" => 400,
    "via1" => 400,
    "metal3" => 400,
    "metal3#e1" => 400,
    "metal3#e2" => 400,
]
thermal_diffisitivities = [
    "core" => 90 / 711 / 2330,
    "box" => 1.38 / 1000 / 2650,
    "clad" => 1.38 / 1000 / 2650,
    "heater" => 400 / 376 / 8960,
    "via2" => 400 / 376 / 8960,
    "metal2" => 400 / 376 / 8960,
    "via1" => 400 / 376 / 8960,
    "metal3" => 400 / 376 / 8960,
    "metal3#e1" => 400 / 376 / 8960,
    "metal3#e2" => 400 / 376 / 8960,
]
# %% [markdown]
# First we load the mesh, get the labels, define a CellField to evaluate the constants on
# and create functions which work on the tags

# %% tags=["remove-stderr", "hide-output"]
model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)
dΩ = Measure(Ω, 1)

labels = get_face_labeling(model)
tags = get_face_tag(labels, num_cell_dims(model))
τ = CellField(tags, Ω)

electrical_conductivity =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in electrical_conductivity)
ϵ_electrical_conductivity(tag) = electrical_conductivity[tag]

thermal_conductivities =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_conductivities)
ϵ_conductivities(tag) = thermal_conductivities[tag]

thermal_diffisitivities =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_diffisitivities)
ϵ_diffisitivities(tag) = thermal_diffisitivities[tag]

# %% [markdown]
# The next step is to define the boundary conditions, this can be done simply via julia-dicts:

# %% tags=["remove-stderr", "hide-output"]
boundary_potentials = Dict(["metal3#e1___None" => 0.4, "metal3#e2___None" => 0.0])
boundary_temperatures = Dict("metal3#e1___None" => 0.4, "metal3#e2___None" => 0.0)

# %% [markdown]
# Now we're ready to do the simulations! First we simulate the electrical potential,
# then we go on with the temperature

# %% tags=["remove-stderr"]
p0 = compute_potential(ϵ_electrical_conductivity ∘ τ, boundary_potentials)
voltage = abs(sum(values(boundary_potentials) .* [-1, 1]))
power = abs(sum(∫(power_density(p0))dΩ))
current = power / voltage
println("Voltage: ", @sprintf("%.2f V", voltage))
println("Current: ", @sprintf("%.2f mA", current * 1e3))
println("Power: ", @sprintf("%.2f mW", power * 1e3))

# %% tags=["remove-stderr"]
T0 = calculate_temperature(ϵ_conductivities ∘ τ, power_density(p0), boundary_temperatures)

# %% [markdown]
# Yay, now we have both fields simulated! Let's get the average temperature of the silicon waveguide.

# %% tags=["remove-stderr"]
Ω_w = Triangulation(model, tags = "core")
dΩ_w = Measure(Ω_w, 1)
silicon_volume = ∑(∫(1)dΩ_w)

println(
    "Average Temperature of silicon waveguide: ",
    ∑(∫(temperature(T0))dΩ_w) / silicon_volume,
    " K",
)

# %% [markdown]
# And we write the fields to a file for visualisation using paraview:

# %% tags=["remove-stderr", "hide-output"]
writevtk(
    Ω,
    "results",
    cellfields = [
        "potential" => potential(p0),
        "current" => current_density(p0),
        "temperature" => temperature(T0),
        "heat flux" => heat_flux(T0),
        "tags" => τ,
    ],
)


# %% [markdown]
# Now let's get to the transient simulation.
# We'll simulate the heatup by starting with zero temperature offset and the power density based on the electrical simulation.
# For the cooldown simulation we start with the steady state temperature profile and no heating.

# %% tags=["remove-stderr"]
figure = Figure()
ax = Axis(
    figure[1, 1],
    ylabel = "Average silicon waveguide temperature / K",
    xlabel = "Time / ms",
)

for (label, power_factor, temperature_factor) in [("heatup", 1, 0), ("cooldown", 0, 1)]
    uₕₜ = calculate_temperature_transient(
        ϵ_conductivities ∘ τ,
        ϵ_diffisitivities ∘ τ,
        power_density(p0) * power_factor,
        boundary_temperatures,
        temperature(T0) * temperature_factor,
        2e-7,
        2e-4,
    )

    #createpvd("poisson_transient_solution_$label") do pvd
    #    for (uₕ, t) in uₕₜ
    #        pvd[t] = createvtk(
    #            Ω,
    #            "poisson_transient_solution_$t" * ".vtu",
    #            cellfields = ["u" => uₕ],
    #        )
    #    end
    #end

    sums = [(t, ∑(∫(u)dΩ_w) / silicon_volume) for (u, t) in uₕₜ]
    t, s = getindex.(sums, 1), getindex.(sums, 2)
    lines!(ax, t * 1e3, s, label = label)
end

axislegend(position = :rc)
display(figure)


# %% [markdown]
# For more complex examples we'd need a iterative solver, i.e. 
#
# ```julia
#   GridapPETSc.with(args = split("-ksp_type cg -pc_type mg -ksp_monitor")) do
#       # Here the code
#   end
# ```
# 
# And then add to each of the calculate-calls the sovler parameter:
#
# ```julia
#   solver = PETScLinearSolver()
# ```
