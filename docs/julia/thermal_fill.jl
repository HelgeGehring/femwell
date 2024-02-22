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
# # Effective thermal conductivity

# Usually, after the design of the circuitry, the empty space on chips gets filled with equally spaced small rectangles in each layer.
# Optically, these structures are supposed to be far enough away that their influence on the structures can be neglected.
# But for thermal considerations, those fill structures can have an impact on the temperature distribution on the chip and thus e.g. on the crosstalk between thermal phase shifters.
# As it's computationally challenging to include all the small cuboids in the model (which is especially for the meshing a major challenge),
# a preferable approach is to consider the filled area as a homogenous area of higher thermal conductivity.
# For this, we calculate the effective thermal conductivity of the filled area by examining a single unit cell.

# To have an intuitively understandable problem, we consider half of the unit cell to be filled with a highly thermally conductive material (metal/silicon) surrounded by a material with low thermal conductance (e.g. silicon dioxide)
# Let's start with the mesh!

# %% 

unitcell_width = 2
unitcell_length = 2
unitcell_thickness = 2

fill_width = 2
fill_length = 1
fill_thickness = 2

# %% tags=["hide-input", "thebe-init", "remove-stderr"]
using PyCall
np = pyimport("numpy")
shapely = pyimport("shapely")
shapely.affinity = pyimport("shapely.affinity")
OrderedDict = pyimport("collections").OrderedDict
meshwell = pyimport("meshwell")
Prism = pyimport("meshwell.prism").Prism
Model = pyimport("meshwell.model").Model

model = Model()

fill_polygon = shapely.box(
    xmin = unitcell_length / 2 - fill_length / 2,
    ymin = unitcell_width / 2 - fill_width / 2,
    xmax = unitcell_length / 2 + fill_length / 2,
    ymax = unitcell_width / 2 + fill_width / 2,
)

fill = Prism(
    polygons = fill_polygon,
    buffers = Dict(0 => 0.0, fill_thickness => 0.0),
    model = model,
    physical_name = "fill",
    mesh_order = 1,
    resolution = Dict("resolution" => 0.2, "SizeMax" => 1.0, "DistMax" => 1.0),
)

unitcell_polygon =
    shapely.box(xmin = 0, ymin = 0, xmax = unitcell_length, ymax = unitcell_width)
unitcell_buffer = Dict(0 => 0.0, unitcell_thickness => 0.0)

unitcell = Prism(
    polygons = unitcell_polygon,
    buffers = unitcell_buffer,
    model = model,
    physical_name = "unitcell",
    mesh_order = 2,
)

"""
BOUNDARIES
"""

right_polygon = shapely.box(
    xmin = unitcell_length,
    ymin = 0,
    xmax = unitcell_length + 1,
    ymax = unitcell_width,
)
right = Prism(
    polygons = right_polygon,
    buffers = unitcell_buffer,
    model = model,
    physical_name = "right",
    mesh_bool = false,
    mesh_order = 0,
)

left_polygon = shapely.box(xmin = -1, ymin = 0, xmax = 0, ymax = unitcell_width)
left = Prism(
    polygons = left_polygon,
    buffers = unitcell_buffer,
    model = model,
    physical_name = "left",
    mesh_bool = false,
    mesh_order = 0,
)

up_polygon = shapely.box(
    xmin = 0,
    ymin = unitcell_width,
    xmax = unitcell_length,
    ymax = unitcell_width + 1,
)
up = Prism(
    polygons = up_polygon,
    buffers = unitcell_buffer,
    model = model,
    physical_name = "up",
    mesh_bool = false,
    mesh_order = 0,
)

down_polygon = shapely.box(xmin = 0, ymin = -1, xmax = unitcell_length, ymax = 0)
down = Prism(
    polygons = down_polygon,
    buffers = unitcell_buffer,
    model = model,
    physical_name = "down",
    mesh_bool = false,
    mesh_order = 0,
)

"""
ASSEMBLE AND NAME ENTITIES
"""
entities = [unitcell, fill, up, down, left, right]

mesh = model.mesh(
    entities_list = entities,
    verbosity = 0,
    filename = "mesh.msh",
    default_characteristic_length = 0.2,
)



# %% tags=["hide-input", "thebe-init", "remove-stderr"]
using Gridap
using GridapGmsh
using Gridap.Geometry
using GridapMakie, CairoMakie

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

# %% [markdown]
# For simplicity, we define the thermal conductivity of the unitcell to be 1 and the thermal conductivity of the fill to be 100, which is almost negligible in comparison.

# %%
thermal_conductivities = ["unitcell" => 1, "fill" => 100.0]

model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)
dΩ = Measure(Ω, 1)

labels = get_face_labeling(model)
tags = get_face_tag(labels, num_cell_dims(model))
τ = CellField(tags, Ω)

thermal_conductivities =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_conductivities)
ϵ_conductivities(tag) = thermal_conductivities[tag]

# %% [markdown]
# We define the temperatures at both sides of the unitcell to define the direction in which we want to estimate the thermal conductivity.
# In other directions the boundaries are considered to be insulating/mirroring.

# %% tags=["remove-stderr", "hide-output"]
boundary_temperatures = Dict("left___unitcell" => 100.0, "right___unitcell" => 0.0)


# %% [markdown]
# We start with calculating the temperature distribution.
# From this, we calculate, analog to how we calculate resistances for electrical simulations first the integral
# $\int σ \left|\frac{\mathrm{d}T}{\mathrm{d}\vec{x}}\right|^2 dA which gives an equivalent of the power
# and calculate from this the effective thermal conductivity.
# We expect the thermal conductivity almost twice as high as the thermal conductivity of the material with lower thermal conductivity.

# %% tags=["remove-stderr"]
T0 = calculate_temperature(
    ϵ_conductivities ∘ τ,
    CellField(x -> 0, Ω),
    boundary_temperatures,
    order = 2,
)
temperature_difference = abs(sum(values(boundary_temperatures) .* [-1, 1]))
power = abs(
    sum(
        ∫(
            (ϵ_conductivities ∘ τ) *
            (norm ∘ gradient(temperature(T0))) *
            (norm ∘ gradient(temperature(T0))),
        )dΩ,
    ),
)

println(
    "Effective thermal conductivity: ",
    (power / temperature_difference^2) /
    ((unitcell_thickness * unitcell_width) / unitcell_length),
)
