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
# # Mode solving anisotropic
# {cite}`Fallahkhair2008`

# %% [markdown]
# ```{caution}
# **This example is under construction**
# As Julia-Dicts are not ordered, the mesh might become incorrect when adjusted (for now, better do the meshing in python)
# ```

# %% tags=["hide-output"]
using PyCall
np = pyimport("numpy")
shapely = pyimport("shapely")
shapely.affinity = pyimport("shapely.affinity")
clip_by_rect = pyimport("shapely.ops").clip_by_rect
OrderedDict = pyimport("collections").OrderedDict
mesh_from_OrderedDict = pyimport("femwell.mesh").mesh_from_OrderedDict

wg_width = 0.8
wg_thickness = 0.6076
core = shapely.geometry.box(-wg_width / 2, 0, +wg_width / 2, wg_thickness)
env = shapely.affinity.scale(core.buffer(5, resolution = 8), xfact = 0.5)

polygons = OrderedDict(
    core = core,
    box = clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad = clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
)

resolutions = Dict("core" => Dict("resolution" => 0.03, "distance" => 0.5))

mesh = mesh_from_OrderedDict(
    polygons,
    resolutions,
    default_resolution_max = 0.1,
    filename = "mesh.msh",
)

# %% tags=["remove-stderr", "hide-output"]
using LinearAlgebra
using SparseArrays

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using Gridap.ReferenceFEs
using GridapGmsh
using GridapMakie, CairoMakie

using Femwell.Maxwell.Waveguide

CairoMakie.inline!(true)

# %% tags=["remove-stderr"]
model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)

labels = get_face_labeling(model)

eye = diagonal_tensor(VectorValue([1.0, 1.0, 1.0]))
n2 = 2.302
Δ = 0.005
yig = TensorValue([n2^2 1im*Δ 0; -1im*Δ n2^2 0; 0 0 n2^2])
epsilons = ["core" => yig, "box" => 1.95^2 * eye, "clad" => 1.0^2 * eye]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]

τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

modes = calculate_modes(model, ε ∘ τ, λ = 1.3, num = 2, order = 1)

plot_mode(modes[1])
plot_mode(modes[2])
modes

# %% [markdown]

# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
