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
# # Mode solving for bent waveguides

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
meshwell = pyimport("meshwell")
meshwell.model = pyimport("meshwell.model")
meshwell.polySurface = pyimport("meshwell.polysurface")

function write_mesh(;
    radius = 2,
    wg_width = 0.5,
    wg_thickness = 0.22,
    sim_left = 0.5,
    sim_right = 4,
    sim_top = 1,
    sim_bottom = 4,
    pml_thickness = 3,
    distance = 1,
)
    core1 = shapely.geometry.box(
        radius - wg_width / 2 - distance / 2,
        0,
        radius + wg_width / 2 - distance / 2,
        wg_thickness,
    )
    core2 = shapely.geometry.box(
        radius - wg_width / 2 + distance / 2,
        0,
        radius + wg_width / 2 + distance / 2,
        wg_thickness,
    )

    env = shapely.geometry.box(
        radius - wg_width / 2 - distance / 2 - sim_left,
        -sim_bottom - pml_thickness,
        radius + wg_width / 2 + distance / 2 + sim_right + pml_thickness,
        wg_thickness + sim_top,
    )

    model = meshwell.model.Model()

    polygons = [
        meshwell.polySurface.PolySurface(
            model = model,
            polygons = core1,
            physical_name = "core1",
            mesh_order = 1,
            resolution = Dict("resolution" => 0.05, "DistMax" => 1.5, "SizeMax" => 0.15),
        ),
        meshwell.polySurface.PolySurface(
            model = model,
            polygons = core2,
            physical_name = "core2",
            mesh_order = 1,
            resolution = Dict("resolution" => 0.05, "DistMax" => 1.5, "SizeMax" => 0.15),
        ),
        meshwell.polySurface.PolySurface(
            model = model,
            polygons = clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
            physical_name = "clad",
            mesh_order = 2,
        ),
        meshwell.polySurface.PolySurface(
            model = model,
            polygons = clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
            physical_name = "box",
            mesh_order = 2,
        ),
    ]

    model.mesh(polygons, filename = "mesh.msh", default_characteristic_length = 0.3)
end

# %% tags=["remove-stderr", "hide-output", "thebe-init"]
using Gridap
using Gridap.Geometry
using Gridap.Visualization
using Gridap.ReferenceFEs
using GridapGmsh
using GridapMakie, CairoMakie

using Femwell.Maxwell.Waveguide

CairoMakie.inline!(true)

# %% tags=["remove-stderr", "hide-output"]
radiuss = 10 # [10.,15.,20.]
wg_width = 0.55
sim_right = 1
sim_bottom = 1
overlaps = Float64[]
Δneff_cylindrical = Float64[]
for radius in radiuss
    write_mesh(
        radius = radius,
        wg_width = wg_width,
        sim_right = sim_right,
        sim_bottom = sim_bottom,
        distance = 0.8,
    )
    model = GmshDiscreteModel("mesh.msh")
    Ω = Triangulation(model)
    labels = get_face_labeling(model)

    τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)
    pml_x = x -> 0.1 * max(0, x[1] - (radius + wg_width / 2 + sim_right))^2
    pml_y = x -> 0.1 * max(0, -x[2] - sim_bottom)^2


    epsilons = ["core1" => 3.48^2, "core2" => 3.48^2, "box" => 1.46^2, "clad" => 1.46^2]
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
    modes_1 = calculate_modes(
        model,
        ε ∘ τ,
        λ = 1.55,
        num = 4,
        order = 1,
        #radius = radius,
        pml = [pml_x, pml_y],
    )
    plot_mode(modes_1[1], absolute = false, field = E)
    plot_mode(modes_1[1], absolute = false, field = H)
    #println(n_eff(modes_1[1]))
    epsilons = ["core1" => 3.48^2, "core2" => 3.48^2, "box" => 1.46^2, "clad" => 1.46^2]
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
    modes = calculate_modes(
        model,
        ε ∘ τ,
        λ = 1.55,
        num = 4,
        order = 1,
        radius = radius,
        pml = [pml_x, pml_y],
    )
    println(n_eff_cylidrical(modes[1]))
    println(n_eff_cylidrical(modes[2]))
    plot_mode(modes[1], absolute = false)
    plot_mode(modes[1], absolute = false, field = H)
    plot_mode(modes[2], absolute = false)
    plot_mode(modes[2], absolute = false, field = H)
    push!(overlaps, abs(overlap(modes_1[1], modes[1])^2))
    push!(Δneff_cylindrical, real(n_eff_cylidrical(modes[1]) - n_eff_cylidrical(modes[2])))

    println("overlaps")
    println(abs(overlap(modes_1[1], modes_1[1])^2))
    println(abs(overlap(modes_1[1], modes_1[2])^2))
    println(abs(overlap(modes_1[1], modes_1[3])^2))
    println(abs(overlap(modes_1[1], modes_1[4])^2))
    println(abs(overlap(modes_1[1], modes[1])^2))
    println(abs(overlap(modes_1[1], modes[2])^2))
    println(abs(overlap(modes[1], modes[1])^2))
    println(abs(overlap(modes[2], modes[2])^2))
    println(abs(overlap(modes[3], modes[3])^2))
    println(abs(overlap(modes[1], modes[2])^2))
    println(abs(overlap(modes[1], modes[3])^2))
    println(abs(overlap(modes[1], modes[4])^2))
end

# %% tags=["hide-output"]

# %% 
# figure = Figure()
# ax = Axis(figure[1, 1], xlabel = "Radius / μm", ylabel = "Overlap single mode->coupled mode / dB")
# lines!(ax, radiuss, 10log10.(overlaps))
# scatter!(ax, radiuss, 10log10.(overlaps), label = "FEM")
# #axislegend()
# display(figure)

# figure = Figure()
# ax = Axis(figure[1, 1], xlabel = "Radius / μm", ylabel = "π/(k_0 * Δneff)")
# lines!(ax, radiuss, π * (2π/1.55 * Δneff_cylindrical).^-1)
# scatter!(ax, radiuss, π * (2π/1.55 * Δneff_cylindrical).^-1, label = "FEM")
# #axislegend()
# display(figure)
