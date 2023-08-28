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

py"""
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import box
from skfem import Basis, ElementDG, ElementTriP0, ElementTriP1
from skfem.io import from_meshio
from tqdm import tqdm

from femwell.culomb import solve_coulomb
from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

core_width = 1.532
electrode_start_x = core_width / 2 + 2.629
electrode_width = 4.4

electrode_left = box(-electrode_start_x - electrode_width, 0.5, -electrode_start_x, 0.5 + 1.8)
electrode_right = box(electrode_start_x, 0.5, electrode_start_x + electrode_width, 0.5 + 1.8)
slab = box(-10, 0, 10, 0.5)
core = box(-core_width / 2, 0.5, core_width / 2, 0.8)
env = slab.buffer(20, resolution=8)

polygons = OrderedDict(
    electrode_left=electrode_left,
    electrode_right=electrode_right,
    slab=slab,
    core=core,
    env=env,
)

resolutions = dict(
    slab={"resolution": 0.1, "distance": 4},
    core={"resolution": 0.1, "distance": 4},
)

mesh = from_meshio(
    mesh_from_OrderedDict(polygons, resolutions, filename="mesh.msh", default_resolution_max=3)
)
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

conductivity = [
    "core" => 7.5,
    "slab" => 28.4,
    "electrode_left" => 3.9,
    "electrode_right" => 3.9,
    "env" => 3.9,
]
σ(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in conductivity)[tag]

p0 = compute_potential(
    σ ∘ τ,
    Dict("electrode_left___slab" => 1.0, "electrode_right___slab" => 0.0),
    order = 2,
)

plot_potential(p0)

fig, _, plt =
    plot(get_triangulation(potential(p0)), (σ ∘ τ) * ∇(potential(p0)) ⋅ VectorValue(1, 0))
Colorbar(fig[1, 2], plt)
display(fig)

voltages = 0:10:100
voltages_neffs = []

for voltage in voltages
    epsilons = [
        "core" => 1.989,
        "env" => 1.445,
        "slab" => 2.211,
        "electrode_left" => 1.445,
        "electrode_right" => 1.445,
    ]
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]

    ε_perturbation =
        0.5 * 2.211^3 * 31e-6 * ∇(potential(p0)) ⋅ VectorValue(-1, 0) *
        voltage *
        ((tag -> get_tag_from_name(labels, "slab") == tag) ∘ τ)

    modes = calculate_modes(
        model,
        (ε ∘ τ + ε_perturbation) * (ε ∘ τ + ε_perturbation),
        λ = 1.55,
        order = 0,
    )
    push!(voltages_neffs, n_eff(modes[1]))
end

# %% tags=["hide-input"]
f = Figure()
ax = Axis(f[1, 1], xlabel = "Voltage / V", ylabel = "Effective refractive index")
lines!(ax, voltages, real(voltages_neffs))
plot!(ax, voltages, real(voltages_neffs))
display(f)

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
