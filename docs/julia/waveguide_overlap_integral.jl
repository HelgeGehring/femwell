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
# # Noise floor of overlap integrals
# 
# The overlap integral between neighboring waveguides is one way to determine the coupling between neighboring waveguides.
# As the modes decay exponentially outside the core after a certain distance, the overlap integral also decays in good approximation exponentially.
# In order to see how good the estimation of the overlap integral is and to see at which distance the coupling is too weak to be estimated precisely using this methodology,
# we calculate the coupling as a function of the distance in the following.
#
# For this example we chose two silicon nitride waveguides with a width of 1μm and a thickness of 0.3μm. Those waveguides are ontop of a silicon dioxide layer and air-clad.

# %% tags=["hide-input","hide-output"]
using PyCall

all_distances = 1:1.0:2
distances = Float64[]
overlaps = Float64[]

mode_1 = nothing
mode_2 = nothing
Ω = nothing

for distance in all_distances
    println(distance)

    py"""
    import numpy as np
    import shapely
    import shapely.geometry
    import shapely.affinity
    from shapely.ops import clip_by_rect, unary_union
    from collections import OrderedDict
    from femwell.mesh import mesh_from_OrderedDict

    distance = $distance
    wg_width = 1
    wg_thickness = 0.3
    core = shapely.geometry.box(-wg_width / 2, 0, +wg_width / 2, wg_thickness)
    core1 = shapely.affinity.translate(core, xoff=distance/2)
    core2 = shapely.affinity.translate(core, xoff=-distance/2)
    env = unary_union([core1, core2]).buffer(max(5,distance), quad_segs = 2)

    polygons = OrderedDict(
        core1=core1,
        core2=core2,
        box = clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
        clad = clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
    )

    resolutions = {f"core{i}": {"resolution": 0.05, "distance": 1} for i in [1,2]}
    mesh = mesh_from_OrderedDict(
        polygons,
        resolutions,
        default_resolution_max = 1,
        filename = "mesh.msh"
    )
    """
    using Gridap
    using Gridap.Geometry
    using Gridap.Visualization
    using Gridap.ReferenceFEs
    using GridapGmsh
    using GridapMakie, CairoMakie

    using Femwell.Maxwell.Waveguide

    order = 1

    CairoMakie.inline!(true)

    model = GmshDiscreteModel("mesh.msh")
    Ω = Triangulation(model)
    labels = get_face_labeling(model)
    τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

    fig = plot(Ω)
    fig.axis.aspect = DataAspect()
    wireframe!(Ω, color = :black, linewidth = 1)
    display(fig)

    epsilons = ["core1" => 1.9963^2, "core2" => 1.444^2, "box" => 1.444^2, "clad" => 1.0^2]
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
    @time modes_1 = calculate_modes(model, ε ∘ τ, λ = 1.55, num = 1, order = order)
    println(n_eff(modes_1[1]))

    if length(distances) == 0
        push!(distances, 0)
        push!(overlaps, overlap(modes_1[1], modes_1[1]))
    end

    epsilons = ["core1" => 1.444^2, "core2" => 1.9963^2, "box" => 1.444^2, "clad" => 1.0^2]
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
    modes_2 = calculate_modes(model, ε ∘ τ, λ = 1.55, num = 1, order = order)
    println(n_eff(modes_2[1]))

    overlap_integral = overlap(modes_1[1], modes_2[1])
    println(overlap_integral)
    push!(distances, distance)
    push!(overlaps, abs(overlap_integral))

    if distance == last(all_distances)
        plot_mode(modes_1[1])
        plot_mode(modes_2[1])
    end
end

# %% tags=["hide-input"]
f = Figure()
ax = Axis(
    f[1, 1],
    title = "Overlap integral between neighboring waveguides",
    xlabel = "Distance [μm]",
    ylabel = "Overlap integral [dB]",
)
lines!(ax, distances, 10 * log10.(overlaps))
plot!(ax, distances, 10 * log10.(overlaps))
display(f)
