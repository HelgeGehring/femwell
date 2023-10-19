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
# # Microstrip Waveguide

# %% [markdown]
# ```{caution}
# **This example is under construction**
# As Julia-Dicts are not ordered, the mesh might become incorrect when adjusted (for now, better do the meshing in python)
# ```

# %% tags=["hide-input"]
using PyCall

microstrip_L = 250e-6
microstrip_w = 10e-6
microstrip_h = 26e-6
microstrip_t = 1e-6
microstrip_σ = 5.8e7
microstrip_ε_r = 11.7

py"""
from collections import OrderedDict

from shapely.geometry import box, LineString

from femwell.mesh import mesh_from_OrderedDict

cell_height = 1.5e-3
cell_width = 1e-3

core = box(-$microstrip_w, $microstrip_h, $microstrip_w, $microstrip_h+$microstrip_t)
si = box(-$microstrip_L, 0, $microstrip_L, $microstrip_h)
air = box(-cell_width, $0, cell_width, cell_height)

bottom = LineString(((-cell_width, 0), (cell_width, 0)))
#bottom_si = LineString(((-$microstrip_L, 1e-6), ($microstrip_L, 1e-6)))
top = LineString(((-cell_width, cell_height), (cell_width, cell_height)))
left = LineString(((-cell_width, 0), (-cell_width, cell_height)))
right = LineString(((cell_width, 0), (cell_width, cell_height)))

polygons = OrderedDict(
    bottom=bottom,
    #bottom_si=bottom_si,
    top=top,
    left=left,
    right=right,
    core=core,
    si=si,
    air=air
)

resolutions = dict(
    core={"resolution": .4e-6, "distance": 40e-6},
    #si={"resolution": 2e-6, "distance": 40e-6},
)

mesh_from_OrderedDict(polygons, resolutions, filename="mesh.msh", default_resolution_max=40e-6)
"""

# %% tags=["remove-stderr"]
using Gridap
using Gridap.Geometry
using GridapGmsh
using Femwell.Maxwell.Waveguide
using Femwell.Maxwell.Electrostatic
using GridapMakie, CairoMakie
import PhysicalConstants.CODATA2018: c_0, μ_0, ε_0
import Unitful: ustrip

model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)
labels = get_face_labeling(model)
τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

f0 = .01e9

epsilons = [
    "core" => 0 - 1im * microstrip_σ / f0 / ustrip(ε_0),
    "si" => microstrip_ε_r,
    "air" => 1,
]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]

modes = calculate_modes(
    model,
    ε ∘ τ,
    λ = 3e8 / f0,
    order = 1,
    metallic_boundaries = ["bottom", "bottom_points", "top", "top_points", "left", "right"],
)

write_mode_to_vtk("mode.vtu", modes[1])
plot_mode(modes[1], vertical = true, same_range = false)
modes

# %% tags=["hide-input"]

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
