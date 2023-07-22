using GLMakie

using Gridap
using GridapGmsh
using Gridap.Geometry

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal


model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)

labels = get_face_labeling(model)
tags = get_face_tag(labels, num_cell_dims(model))
τ = CellField(tags, Ω)

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

boundary_conditions = Dict(["metal3#e1___None" => 0.1, "metal3#e2___None" => 0.0])

p0 = compute_potential(ϵ_electrical_conductivity ∘ τ, boundary_conditions)

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

T0 = calculate_temperature(ϵ_conductivities ∘ τ, power_density(p0), temperatures)

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
    power_density(p0) * 0,
    temperatures,
    temperature(T0),
    .5e-6,
    2e-5,
)

createpvd("poisson_transient_solution") do pvd
    for (uₕ, t) in uₕₜ
        pvd[t] =
            createvtk(Ω, "poisson_transient_solution_$t" * ".vtu", cellfields = ["u" => uₕ])
    end
end

Ω_w = Triangulation(model, tags = "core")
dΩ_w = Measure(Ω_w, 1)
sums = [(t, ∑(∫(u)dΩ_w)) for (u, t) in uₕₜ]
lines(sums)
