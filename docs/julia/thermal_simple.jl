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
# # Mode solving

# %% tags=["hide-input", "thebe-init"]
using Gridap
using Gridap.Geometry

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

# %% tags=["hide-output"]
domain = (0, 1, 0, 1)
partition = (20, 20)
model = CartesianDiscreteModel(domain, partition)
labels = get_face_labeling(model)
add_tag!(labels, "left", [1, 2, 5])
add_tag!(labels, "right", [3, 4, 6])
tags = get_face_tag(labels, num_cell_dims(model))
Ω = Triangulation(model)
dΩ = Measure(Ω, 1)
τ = CellField(tags, Ω)
constant = tag -> 1

# %% tags=["hide-input"]
p0 = compute_potential(constant ∘ τ, Dict("left" => 0.0, "right" => 1.0))

average_potential = ∑(∫(potential(p0))dΩ) / ∑(∫(1)dΩ)
println(average_potential)

# %% tags=["hide-input"]
T0 = calculate_temperature(constant ∘ τ, power_density(p0), Dict("boundary" => 0.0))

writevtk(
    Ω,
    "results",
    cellfields = [
        "potential" => potential(p0),
        "current" => current_density(p0),
        "temperature" => temperature(T0),
    ],
)
# Dict{String,Float64}()

# %% tags=["hide-input"]
T_transient = calculate_temperature_transient(
    constant ∘ τ,
    constant ∘ τ,
    power_density(p0) * 0,
    Dict("boundary" => 0.0),
    temperature(T0),
    1e-4,
    1.e-3,
)
sums = [(t, ∑(∫(u)dΩ) / ∑(∫(1)dΩ)) for (u, t) in T_transient]
println(sums)

# %% tags=["hide-input"]
T_transient = calculate_temperature_transient(
    constant ∘ τ,
    constant ∘ τ,
    power_density(p0),
    Dict{String,Float64}(),
    temperature(T0) * 0,
    1e-1,
    1.0,
)
sums = [(t, ∑(∫(u)dΩ) / ∑(∫(1)dΩ)) for (u, t) in T_transient]
println(sums)
