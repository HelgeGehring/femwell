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
# # 3D thermal phase shifter

# %% tags=["remove-stderr", "hide-input", "hide-output"]
using GridapMakie
using Printf

include("./js.jl")
using .JS: FigureContainer

using Gridap
using GridapGmsh
using Gridap.Geometry
using GridapPETSc

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

using IJulia
using CairoMakie
using Markdown

using WGLMakie
WGLMakie.activate!()

Makie.inline!(true) # Make sure to inline plots into Documenter output!

dir = @__DIR__
# run(`python $dir/heater_3d_mesh.py`)

# %% tags=["remove-stderr", "hide-input", "hide-output"]

function visualize_mesh(model)
    volume_tags = [
        "box",
        "core",
        "heater",
        "via2",
        "via1",
        "metal3",
        "metal2",
        "metal3#e1",
        "metal3#e2",
    ]
    colors = [:grey, :black, :orange, :green, :green, :brown, :brown, :blue, :blue]
    alphas = [1.0, 1.0, 0.2, 0.5, 0.5, 0.1, 0.1, 0.75, 0.75]
    fig = Figure(fontsize = 16)
    xlims = (-20e-6, 60e-6)
    ylims = (-5e-6, 5e-6)
    zlims = (-5e-6, 5e-6)
    aspect = (xlims[2] - xlims[1]) / (ylims[2] - ylims[1])
    tickfunc = values -> (x -> string(x * 1e6)).(values)

    # Mesh
    ax =
        fig[1, 2] = Axis3(
            fig[1, 2],
            aspect = (aspect, 1.0, 1.0),
            xlabel = "x (μm)",
            ylabel = "\ny (μm)",
            zlabel = "z (μm)\n\n",
            title = "Heater Mesh",
            limits = (xlims, ylims, zlims),
            titlesize = 28,
            xtickformat = tickfunc,
            ytickformat = tickfunc,
            ztickformat = tickfunc,
            azimuth = 51 / 40 * pi + pi / 20,
            elevation = pi / 8 + pi / 10,
        )
    Camera3D(ax.scene)

    for (i, tag) in enumerate(volume_tags)
        ∂Ω_region = BoundaryTriangulation(model, tags = [tag])
        wireframe!(
            fig[1, 2],
            ∂Ω_region,
            color = colors[i],
            alpha = alphas[i],
            linewidth = 0.1,
            fontsize = 16,
        )
    end
    wf = [PolyElement(color = color) for color in colors]
    Legend(
        fig[1, 1],
        wf,
        volume_tags,
        halign = :left,
        valign = :top,
        framevisible = false,
        labelsize = 16,
    )

    display(FigureContainer(fig))
end

function visualize_field(model, field, label; plotargs...)
    volume_tags = [
        "box",
        "core",
        "heater",
        "via2",
        "via1",
        "metal3",
        "metal2",
        "metal3#e1",
        "metal3#e2",
    ]
    fontsize = 10
    fig = Figure(fontsize = fontsize, resolution = (800, 500))
    xlims = (-20e-6, 60e-6)
    ylims = (-5e-6, 5e-6)
    zlims = (-5e-6, 5e-6)
    aspect = (xlims[2] - xlims[1]) / (ylims[2] - ylims[1])
    tickfunc = values -> (x -> string(x * 1e6)).(values)

    # Setup plot
    Ω = Triangulation(model, tags = volume_tags)
    ax =
        fig[1, 1] = Axis3(
            fig[1, 1],
            aspect = (aspect, 1.0, 1.0),
            xlabel = "x (μm)",
            ylabel = "y (μm)",
            zlabel = "z (μm)",
            title = label,
            limits = (xlims, ylims, zlims),
            titlesize = 18,
            xtickformat = tickfunc,
            ytickformat = tickfunc,
            ztickformat = tickfunc,
            azimuth = 51 / 40 * pi + pi / 20,
            elevation = pi / 8 + pi / 10,
        )
    Camera3D(ax.scene)

    # Plot
    plt = plot!(fig[1, 1], Ω, field, shading = true; plotargs...)
    Colorbar(fig[1, 2], plt, label = label, vertical = true, flipaxis = false)

    display(FigureContainer(fig))
end

function plot_time(uₕₜ::Gridap.ODEs.TransientFETools.TransientFESolution; label = "")
    # Set up figure
    volume_tags = [
        "box",
        "core",
        "heater",
        "via2",
        "via1",
        "metal3",
        "metal2",
        "metal3#e1",
        "metal3#e2",
    ]
    fig = CairoMakie.Figure()
    xlims = (-20e-6, 60e-6)
    ylims = (-5e-6, 5e-6)
    zlims = (-5e-6, 5e-6)
    aspect = (xlims[2] - xlims[1]) / (ylims[2] - ylims[1])
    tickfunc = values -> (x -> string(x * 1e6)).(values)
    Ω = Triangulation(model, tags = volume_tags)

    # Set up axis
    it = iterate(uₕₜ)
    uh = Any[nothing]
    time = Any[nothing]
    state = Any[nothing]
    (uh[1], time[1]), state[1] = it
    t_uh_time = Observable(time[1])
    timetitle = lift(time -> label * @sprintf("   t = %.2f ms", time * 1e3), t_uh_time)
    Axis3(
        fig[1, 1],
        aspect = (aspect, 1.0, 1.0),
        xlabel = "x (μm)",
        ylabel = "\ny (μm)",
        zlabel = "z (μm)",
        title = timetitle,
        limits = (xlims, ylims, zlims),
        xtickformat = tickfunc,
        ytickformat = tickfunc,
        ztickformat = tickfunc,
        azimuth = 51 / 40 * pi + pi / 20,
        elevation = pi / 8 + pi / 10,
    )

    # Plot new iteration
    u = lift(t_uh_time) do t
        while (time[1] < t) || (time[1] ≈ t)
            if time[1] ≈ t
                return uh[1]
            end
            it = iterate(uₕₜ, state[1])
            if isnothing(it)
                return uh[1]
            end
            (uh[1], time[1]), state[1] = it
        end
        return uh[1]
    end

    # Create the actual plot
    plt = plot!(fig[1, 1], Ω, u, colorrange = (0, 110.0), colormap = Reverse(:RdBu))
    Colorbar(fig[1, 2], plt, vertical = true, label = "Temperature (K)")
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    t_uh_time, fig
end

# %% tags=["remove-stderr"]
model = GmshDiscreteModel("$dir/mesh.msh")
visualize_mesh(model)

# %% [markdown]
# Let's start with defining the constants for the simulation:

# %% tags=["hide-output"]
electrical_conductivity = [
    "core" => 1e-6,
    "box" => 1e-13,
    "clad" => 1e-13,
    "heater" => 2.3e6,
    "via2" => 3.77e7,
    "metal2" => 3.77e7,
    "via1" => 3.77e7,
    "metal3" => 3.77e7,
    "metal3#e1" => 3.77e7,
    "metal3#e2" => 3.77e7,
]
thermal_conductivities = [
    "core" => 90,
    "box" => 1.38,
    "clad" => 1.38,
    "heater" => 400,
    "via2" => 400,
    "metal2" => 400,
    "via1" => 400,
    "metal3" => 400,
    "metal3#e1" => 400,
    "metal3#e2" => 400,
]
thermal_diffisitivities = [
    "core" => 90 / 711 / 2330,
    "box" => 1.38 / 1000 / 2650,
    "clad" => 1.38 / 1000 / 2650,
    "heater" => 400 / 376 / 8960,
    "via2" => 400 / 376 / 8960,
    "metal2" => 400 / 376 / 8960,
    "via1" => 400 / 376 / 8960,
    "metal3" => 400 / 376 / 8960,
    "metal3#e1" => 400 / 376 / 8960,
    "metal3#e2" => 400 / 376 / 8960,
]
# %% [markdown]
# First we load the mesh, get the labels, define a CellField to evaluate the constants on
# and create functions which work on the tags

# %% tags=["remove-stderr", "hide-output"]
Ω = Triangulation(model)
dΩ = Measure(Ω, 1)

labels = get_face_labeling(model)
tags = get_face_tag(labels, num_cell_dims(model))
τ = CellField(tags, Ω)

electrical_conductivity =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in electrical_conductivity)
ϵ_electrical_conductivity(tag) = electrical_conductivity[tag]

thermal_conductivities =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_conductivities)
ϵ_conductivities(tag) = thermal_conductivities[tag]

thermal_diffisitivities =
    Dict(get_tag_from_name(labels, u) => v for (u, v) in thermal_diffisitivities)
ϵ_diffisitivities(tag) = thermal_diffisitivities[tag]

# %% [markdown]
# The next step is to define the boundary conditions, this can be done simply via julia-dicts:

# %% tags=["remove-stderr", "hide-output"]
boundary_potentials = Dict(["metal3#e1___None" => 0.4, "metal3#e2___None" => 0.0])
boundary_temperatures =
    Dict("metal3#e1___None" => 0.4, "metal3#e2___None" => 0.0, "box___None" => 0.0)

# %% [markdown]
# Now we're ready to do the simulations! First we simulate the electrical potential,
# then we go on with the temperature

# %% tags=["remove-stderr"]
p0 = compute_potential(ϵ_electrical_conductivity ∘ τ, boundary_potentials)
voltage = abs(sum(values(boundary_potentials) .* [-1, 1]))
power = abs(sum(∫(power_density(p0))dΩ))
current = power / voltage
println("Voltage: ", @sprintf("%.2f V", voltage))
println("Current: ", @sprintf("%.2f mA", current * 1e3))
println("Power: ", @sprintf("%.2f mW", power * 1e3))

# %% tags=["remove-stderr"]
T0 = calculate_temperature(ϵ_conductivities ∘ τ, power_density(p0), boundary_temperatures)

# %% [markdown]
# Yay, now we have both fields simulated! Let's get the average temperature of the silicon waveguide.

# %% tags=["remove-stderr"]
Ω_w = Triangulation(model, tags = "core")
dΩ_w = Measure(Ω_w, 1)
silicon_volume = ∑(∫(1)dΩ_w)

println(
    "Average Temperature of silicon waveguide: ",
    ∑(∫(temperature(T0))dΩ_w) / silicon_volume,
    " K",
)

# %% tags=["remove-stderr"]
visualize_field(
    model,
    potential(p0),
    "Potential";
    colormap = :viridis,
    colorrange = (0, 0.5),
    fontsize = 14,
)

# %% tags=["remove-stderr"]
visualize_field(
    model,
    current_density(p0),
    "Current Density";
    colormap = :viridis,
    colorrange = (0, 1.75e10),
    fontsize = 14,
)

# %% tags=["remove-stderr"]
visualize_field(
    model,
    temperature(T0),
    "Temperature";
    colormap = Reverse(:RdBu),
    colorrange = (0, 120.0),
    fontsize = 14,
)

# %% tags=["remove-stderr"]
visualize_field(
    model,
    heat_flux(T0),
    "Heat Flux";
    colormap = Reverse(:RdBu),
    colorrange = (3.0e4, 4e9),
    fontsize = 14,
)

# %% [markdown]
# Alternatively, we write the fields to a file for visualisation using paraview:

# %% tags=["remove-stderr", "hide-output"]
writevtk(
    Ω,
    "$dir/results",
    cellfields = [
        "potential" => potential(p0),
        "current" => current_density(p0),
        "temperature" => temperature(T0),
        "heat flux" => heat_flux(T0),
        "tags" => τ,
    ],
)


# %% [markdown]
# Now let's get to the transient simulation.
# We'll simulate the heatup by starting with zero temperature offset and the power density based on the electrical simulation.
# For the cooldown simulation we start with the steady state temperature profile and no heating.

# %% tags=["remove-stderr"]
CairoMakie.activate!()
figure = Figure()
ax = Axis(
    figure[1, 1],
    ylabel = "Average silicon waveguide temperature / K",
    xlabel = "Time / ms",
)

# %% [markdown]
# Here we set the parameters for the transient simulation
# %% tags=["remove-stderr", "hide-output"]
Δt = 2e-7
t_end = 2e-4
plot_every = 200
timestamps = collect(Δt:Δt*plot_every:t_end)
total_render_time = 5.0
fps = ceil(Int, length(timestamps) / total_render_time)
html_sources = String[]
lins = CairoMakie.Lines[]
for (label, power_factor, temperature_factor) in [("heatup", 1, 0), ("cooldown", 0, 1)]
    uₕₜ = calculate_temperature_transient(
        ϵ_conductivities ∘ τ,
        ϵ_diffisitivities ∘ τ,
        power_density(p0) * power_factor,
        boundary_temperatures,
        temperature(T0) * temperature_factor,
        Δt,
        t_end,
    )

    # Create animation
    t_uh_time, fig = plot_time(uₕₜ; label)
    record(fig, "$dir/animation-$label.gif", timestamps, framerate = fps) do this_t
        t_uh_time[] = this_t
    end

    # Optionally write the solution to a file
    #createpvd("poisson_transient_solution_$label") do pvd
    #    for (uₕ, t) in uₕₜ
    #        pvd[t] = createvtk(
    #            Ω,
    #            "poisson_transient_solution_$t" * ".vtu",
    #            cellfields = ["u" => uₕ],
    #        )
    #    end
    #end

    sums = [(t, ∑(∫(u)dΩ_w) / silicon_volume) for (u, t) in uₕₜ]
    t, s = getindex.(sums, 1), getindex.(sums, 2)
    push!(lins, lines!(ax, t * 1e3, s, label = label))
end

# %% [markdown]
# If we display the animations we can visualize the heatup and cooldown of the device:
# ![heatup animation](animation-heatup.gif)
# ![cooldown animation](animation-cooldown.gif)

# %% [markdown]
# Finally, we can plot an average of the temperature over time: 
# %% tags=["remove-stderr"]
Legend(
    figure[1, 1],
    lins,
    ["heatup", "cooldown"],
    halign = :right,
    valign = :center,
    labelsize = 16,
)
display(figure)


# %% [markdown]
# For more complex examples we'd need a iterative solver, i.e. 
#
# ```julia
#   GridapPETSc.with(args = split("-ksp_type cg -pc_type mg -ksp_monitor")) do
#       # Here the code
#   end
# ```
# 
# And then add to each of the calculate-calls the solver parameter:
#
# ```julia
#   solver = PETScLinearSolver()
# ```
