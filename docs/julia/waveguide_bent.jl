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
OrderedDict = pyimport("collections").OrderedDict
mesh_from_OrderedDict = pyimport("femwell.mesh").mesh_from_OrderedDict



function write_mesh(;
    radius = 2,
    wg_width = 0.5,
    wg_thickness = 0.22,
    sim_left = 0.5,
    sim_right = 4,
    sim_top = 1,
    sim_bottom = 4,
    pml_thickness = 3,
)
    core =
        shapely.geometry.box(radius - wg_width / 2, 0, radius + wg_width / 2, wg_thickness)

    env = shapely.geometry.box(
        radius - wg_width / 2 - sim_left,
        -sim_bottom - pml_thickness,
        radius + wg_width / 2 + sim_right + pml_thickness,
        wg_thickness + sim_top,
    )

    polygons = OrderedDict(
        core = core,
        box = clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
        clad = clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
    )

    resolutions = Dict(
        "core" => Dict("resolution" => 0.01, "distance" => 1.5),
        "slab" => Dict("resolution" => 0.01, "distance" => 1.5),
    )

    mesh_from_OrderedDict(
        polygons,
        resolutions,
        default_resolution_max = 0.2,
        filename = "mesh.msh",
    )
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

# %% tags=["remove-stderr"]
radiuss = 1:0.5:5
neffs = ComplexF64[]
for radius in radiuss
    write_mesh(radius = radius)
    model = GmshDiscreteModel("mesh.msh")
    Ω = Triangulation(model)
    labels = get_face_labeling(model)

    epsilons = ["core" => 3.48^2, "box" => 1.46^2, "clad" => 1.0^2]
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]

    τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)
    pml_x = x -> 0.1 * max(0, x[1] - (radius + wg_width / 2 + sim_right))^2
    pml_y = x -> 0.1 * max(0, -x[2] - sim_bottom)^2

    modes = calculate_modes(model, ε ∘ τ, λ = 1.55, order = 1)
    println(n_eff(modes[1]))
    modes = calculate_modes(
        model,
        ε ∘ τ,
        λ = 1.55,
        order = 1,
        radius = radius,
        pml = [pml_x, pml_y],
        k0_guess = modes[1].k,
    )
    println(n_eff(modes[1]))
    println(log10(abs(imag(n_eff(modes[1])))))
    plot_mode(modes[1], absolute = true)
    push!(neffs, n_eff(modes[1]))
end

display(neffs)

# %% 
radiuss_reference = 1:5
neff_fd = [2.40762, 2.39421, 2.39204, 2.39128, 2.39091]
neff_fmm = [2.40791, 2.39433, 2.39224, 2.39142, 2.39093]
log10imags_fmm = [-3.63721, -6.48982, -9.30488, -9.97048, -10.36615]

# %% 
figure = Figure()
ax = Axis(figure[1, 1], xlabel = "Radius / μm", ylabel = "log(-imag(n_eff))")
lines!(ax, radiuss, log10.(-imag(neffs)))
scatter!(ax, radiuss, log10.(-imag(neffs)), label = "FEM")
scatter!(ax, radiuss_reference, log10imags_fd, label = "FD")
scatter!(
    ax,
    radiuss_reference,
    log10imags_fmm,
    label = "FMM",
    color = :red,
    marker = :xcross,
)
axislegend()
display(figure)

# %% 
figure = Figure()
ax = Axis(figure[1, 1], xlabel = "Radius / μm", ylabel = "real(n_eff)")
lines!(ax, radiuss, real(neffs))
scatter!(ax, radiuss, real(neffs), label = "FEM")
scatter!(ax, radiuss_reference, log10imags_fd, label = "FD")
scatter!(
    ax,
    radiuss_reference,
    log10imags_fmm,
    label = "FMM",
    color = :red,
    marker = :xcross,
)
axislegend()
display(figure)
