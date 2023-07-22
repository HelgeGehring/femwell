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
# # Coupled mode theory

# %% [markdown]
# ```{caution}
# **This example is under construction**
# ```

# %% tags=["hide-output"]
using PyCall
np = pyimport("numpy")
shapely = pyimport("shapely")
shapely.affinity = pyimport("shapely.affinity")
clip_by_rect = pyimport("shapely.ops").clip_by_rect
OrderedDict = pyimport("collections").OrderedDict
mesh_from_OrderedDict = pyimport("femwell.mesh").mesh_from_OrderedDict

wg_width = 2.5
wg_thickness = 0.3
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
    default_resolution_max = 10,
    filename = "mesh.msh",
)

# %%
using Gridap
using Gridap.Geometry
using Gridap.Visualization
using Gridap.ReferenceFEs
using GridapGmsh
using GridapMakie, GLMakie

using Femwell.Maxwell.Waveguide

GLMakie.inline!(true)

model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)
#fig = plot(Ω)
#fig.axis.aspect=DataAspect()
#wireframe!(Ω, color=:black, linewidth=1)
#display(fig)

labels = get_face_labeling(model)

epsilons = ["core" => 1.9963^2, "box" => 1.444^2, "clad" => 1.0^2]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]


#dΩ = Measure(Ω, 1)
τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

modes = calculate_modes(ε ∘ τ, λ = 1.55, num = 2, order = 1)
println(n_eff(modes[1]))
write_mode_to_vtk("mode", modes[1])
modes

#epsilons_p = ["core"=>-.1im,"box"=>-.0im,"clad"=>-.0im]
#ε_p(tag) = Dict(get_tag_from_name(labels, u)=>v for (u,v) in epsilons_p)[tag]
#println(perturbed_neff(modes[1],ε_p∘τ))
#println((perturbed_neff(modes[1],ε_p∘τ)-n_eff(modes[1]))*2+n_eff(modes[1]))

#modes = calculate_modes(ε∘τ + ε_p∘τ, λ=1.55, num=6, radius=3, order=2)
#println(n_eff(modes[1]))


#plot_mode(modes[1])
