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
# # Mode solving

# %% [markdown]
# ```{caution}
# **This example is under construction**
# As Julia-Dicts are not ordered, the mesh might become incorrect when adjusted (for now, better do the meshing in python)
# ```

# %% tags=["hide-output", "thebe-init"]
using PyCall
np = pyimport("numpy")
shapely = pyimport("shapely")
shapely.affinity = pyimport("shapely.affinity")
clip_by_rect = pyimport("shapely.ops").clip_by_rect
OrderedDict = pyimport("collections").OrderedDict
mesh_from_OrderedDict = pyimport("femwell.mesh").mesh_from_OrderedDict

wg_width = 0.5
wg_thickness = 0.22
slab_width = 3
slab_thickness = 0.11
core = shapely.geometry.box(-wg_width / 2, 0, wg_width / 2, wg_thickness)
slab = shapely.geometry.box(slab_width / 2, 0, slab_width / 2, slab_thickness)
env = shapely.affinity.scale(core.buffer(3, resolution = 8), xfact = 0.5)

polygons = OrderedDict(
    core = core,
    slab = slab,
    box = clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad = clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
)

resolutions = Dict(
    "core" => Dict("resolution" => 0.03, "distance" => 0.5),
    "slab" => Dict("resolution" => 0.03, "distance" => 0.5),
)

mesh = mesh_from_OrderedDict(
    polygons,
    resolutions,
    default_resolution_max = 1,
    filename = "mesh.msh",
)

# %% tags=["remove-stderr", "hide-output", "thebe-init"]
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
#fig = plot(Ω)
#fig.axis.aspect=DataAspect()
#wireframe!(Ω, color=:black, linewidth=1)
#display(fig)

labels = get_face_labeling(model)

epsilons = ["core" => 3.5^2, "slab" => 1.44^2, "box" => 1.444^2, "clad" => 1.44^2]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]


#dΩ = Measure(Ω, 1)
τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

modes = calculate_modes(model, ε ∘ τ, λ = 1.55, num = 1, order = 1)
println(n_eff(modes[1]))
# write_mode_to_vtk("mode", modes[2])

plot_mode(modes[1])
#plot_mode(modes[2])
modes

# %% [markdown]
# ## Perturbations
# Let's add a minor perturbation and calculate the effective refractive index using perturbation theory:

# %% tags=["remove-stderr"]
epsilons_p = ["core" => -0.1im, "box" => -0.0im, "slab" => -0.0im, "clad" => -0.0im]
ε_p(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons_p)[tag]

println(perturbed_neff(modes[1], ε_p ∘ τ))

# %% [markdown]
# For comparison, we also calculate directly the effective refractive index:

# %% tags=["remove-stderr"]

modes_p = calculate_modes(model, ε ∘ τ + ε_p ∘ τ, λ = 1.55, num = 2, order = 1)
println(n_eff(modes_p[1]))
plot_mode(modes_p[1])

# %%
