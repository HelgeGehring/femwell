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

unitcell_width = 2
unitcell_length = 2
unitcell_thickness = 2

fill_width = 2
fill_length = 1
fill_thickness = 2

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


# %%
thermal_conductivities = ["unitcell" => 1, "fill" => 100.0]
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

thermal_conductivities =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_conductivities)
ϵ_conductivities(tag) = thermal_conductivities[tag]

# %% [markdown]
# The next step is to define the boundary conditions, this can be done simply via julia-dicts:

# %% tags=["remove-stderr", "hide-output"]
boundary_temperatures = Dict("left___unitcell" => 100.0, "right___unitcell" => 0.0)

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
#println(∑(∫(1)dΩ))
(power / temperature_difference^2) /
((unitcell_thickness * unitcell_width) / unitcell_length)
