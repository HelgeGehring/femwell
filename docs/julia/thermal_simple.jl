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
# # Match theoretical model for electro-optic simulation

# %% tags=["hide-input", "thebe-init"]
using Gridap
using Gridap.Geometry
using GridapMakie, CairoMakie

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

# %% [markdown]
# We start with setting up a square domain.
# For the boundary conditions, we tag the left and the right side of the model.
# Furthermore, we create a function which returns 1 indipendent of the tag which is the parameter to descrie the constants of the simplified model.

# %% tags=["hide-output", "remove-stderr"]
domain = (0, 1.0, 0, 1.0)
partition = (20, 20)
model = simplexify(CartesianDiscreteModel(domain, partition))
labels = get_face_labeling(model)
add_tag!(labels, "left", [7])
add_tag!(labels, "right", [8])
tags = get_face_tag(labels, num_cell_dims(model))
Ω = Triangulation(model)
dΩ = Measure(Ω, 1)
τ = CellField(tags, Ω)
constant = tag -> 1

# %% [markdown]
# ## Electrostatic
# The first step ist to calculate the potential (assuming the electrical resistivity / conductivity to be ρ=1).
# For this we solve the electrostatic equation $Δϕ = 0$ and define the voltage at two oppositing boundaries to 0V at $x=0$ and 1V at $x=1$.
# The theoretical solution of this function is a linear function.
# $$
#   ϕ(x)=x
# $$
# This would mean the average of the potential over the domain should be
# $$
#   \int ϕ dA / \int 1 dA = 0.5
# $$

# %% tags=[]
p0 = compute_potential(constant ∘ τ, Dict("left" => 0.0, "right" => 1.0))
fig, _, plt = plot(Ω, potential(p0), colormap = :cool)
Colorbar(fig[1, 2], plt)
display(fig)

# %% tags=[]
average_potential = ∑(∫(potential(p0))dΩ) / ∑(∫(1)dΩ)
println("The computed value for the average potential is $average_potential")

# %% [markdown]
# The current density can be calculated as
# $$
# i = ρ \frac{\mathrm{d}ϕ}{\mathrm{d}ϕ} = 1 
# $$
# and thus the averaged current density over the domain to be also 1.

average_current_density = ∑(∫(current_density(p0))dΩ) / ∑(∫(1)dΩ)
println("The computed value for the average current density is $average_current_density")

# %% tags=["hide-input"]
T0 = calculate_temperature(constant ∘ τ, power_density(p0), Dict("boundary" => 0.0))

# %% tags=["hide-input"]
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
    1e-3,
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
