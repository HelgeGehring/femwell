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
# # Lithium niobate phase-shifter

# Reproducing {cite}`Han2022`

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

core = box(-$microstrip_w, $microstrip_h, $microstrip_w, $microstrip_h+$microstrip_t)
si = box(-$microstrip_L, 0, $microstrip_L, $microstrip_h)
air = box(-$microstrip_L, $microstrip_h, $microstrip_L, $microstrip_h*3)

bottom = LineString(((-$microstrip_L, 0), ($microstrip_L, 0)))

polygons = OrderedDict(
    bottom=bottom,
    core=core,
    si=si,
    air=air
)

resolutions = dict(
    core={"resolution": .1e-6, "distance": 4e-6},
)

mesh_from_OrderedDict(polygons, resolutions, filename="mesh.msh", default_resolution_max=5e-6)
"""

# %% tags=["remove-stderr"]
using Gridap
using Gridap.Geometry
using GridapGmsh
using Femwell.Maxwell.Waveguide
using Femwell.Maxwell.Electrostatic
using GridapMakie, CairoMakie

model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)
labels = get_face_labeling(model)
τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)


epsilons = ["core" => 1 + 1im * microstrip_σ / 1e9, "si" => microstrip_ε_r, "air" => 1]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]

modes = calculate_modes(
    model,
    ε ∘ τ,
    λ = 3e8 / 1e9,
    order = 1,
    metallic_boundaries = ["bottom"],
)

plot_mode(modes[1])
modes

# %% tags=["hide-input"]

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
