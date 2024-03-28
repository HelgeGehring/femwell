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
#     name: julia-1.10
# ---

# %% [markdown]
# # Match theoretical model for electro-optic simulation

# %% tags=["hide-input", "thebe-init", "remove-stderr"]
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
domain = (-1.0, 1.0, -1.0, 1.0)
partition = (20, 20)
model = simplexify(CartesianDiscreteModel(domain, partition))
labels = get_face_labeling(model)
add_tag!(labels, "left", [1, 3, 7])
add_tag!(labels, "right", [2, 4, 8])
tags = get_face_tag(labels, num_cell_dims(model))
Ω = Triangulation(model)
dΩ = Measure(Ω, 1)
τ = CellField(tags, Ω)
constant_21 = tag -> 21
constant_42 = tag -> 42

# %% [markdown]
# ## Electrostatic
# The first step ist to calculate the potential (assuming the electrical conductivity to be k=42).
# For this we solve the electrostatic equation $Δϕ = 0$ and define the voltage at two oppositing boundaries to 0V at $x=0$ and 1V at $x=1$.
# The theoretical solution of this function is a linear function.
#
# $$
#   ϕ(x)=x
# $$
#
# This would mean the average of the potential over the domain should be
#
# $$
#   \int ϕ dA / \int 1 dA = 0
# $$

# %% tags=[]
p0 = compute_potential(constant_42 ∘ τ, Dict("left" => -1.0, "right" => 1.0))
fig, _, plt = plot(Ω, potential(p0), colormap = :cool)
Colorbar(fig[1, 2], plt)
display(fig)

# %% tags=[]
average_potential = ∑(∫(potential(p0))dΩ) / ∑(∫(1)dΩ)
println("The computed value for the average potential is $average_potential")

# %% [markdown]
# The current density can be calculated as
#
# $$
# i = k \frac{\mathrm{d}ϕ}{\mathrm{d}ϕ} = 42
# $$
#
# and thus the averaged current density over the domain is 42.

average_current_density = ∑(∫(current_density(p0))dΩ) / ∑(∫(1)dΩ)
println("The computed value for the average current density is $average_current_density")


# %% [markdown]
# Using this value, we can caluclate the average power density as
#
# $$
# p = k i^2
# $$
#
# and thus the averaged power density over the domain is also 42.

# %% tags=[]
average_power_density = ∑(∫(power_density(p0))dΩ) / ∑(∫(1)dΩ)
println("The computed value for the average current density is $average_power_density")

# %% [markdown]
# ## Thermal steady state simulation
# Now we calculate the thermal steady state based on the previously calculated locally applied power.
# For this we chose the thermal conductivity to be $k_{thermal}=21$ and set the boundaries to 0.

# %% tags=[]
T0 = calculate_temperature(constant_21 ∘ τ, power_density(p0), Dict("boundary" => 0.0))
T_average = ∑(∫(temperature(T0))dΩ) / ∑(∫(1)dΩ)
fig, _, plt = plot(Ω, temperature(T0), colormap = :hot)
Colorbar(fig[1, 2], plt)
display(fig)

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

# %% [markdown]
# ## Thermal transient state simulation
# For the simulation of the transient starting with the steady state solution we expect the temperature not to change.
# Also, we don't expect it to depend on the thermal thermal diffusitivity.

# %% tags=[]
T_transient = calculate_temperature_transient(
    constant_21 ∘ τ,
    constant_42 ∘ τ,
    power_density(p0),
    Dict("boundary" => 0.0),
    temperature(T0),
    1e-5,
    1e-3,
)
sums = [(t, ∑(∫(u)dΩ) / ∑(∫(1)dΩ) / T_average) for (u, t) in T_transient]
display(lines(sums))

# %% tags=[]
T_transient = calculate_temperature_transient(
    constant_21 ∘ τ,
    constant_42 ∘ τ,
    power_density(p0) * 0,
    Dict("boundary" => 0.0),
    temperature(T0),
    1e-5,
    2e-2,
)
sums = [(t, ∑(∫(u)dΩ) / ∑(∫(1)dΩ) / T_average) for (u, t) in T_transient]
display(lines(sums))

# %% tags=[]
T_transient = calculate_temperature_transient(
    constant_21 ∘ τ,
    constant_42 ∘ τ,
    power_density(p0),
    Dict{String,Float64}(),
    temperature(T0) * 0,
    1e-2,
    1.0,
)
sums = [(t, ∑(∫(u)dΩ) / ∑(∫(1)dΩ)) for (u, t) in T_transient]
display(lines(sums))
