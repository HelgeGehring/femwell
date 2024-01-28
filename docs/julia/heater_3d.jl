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
# # 3D thermal phase shifter

# %% tags=["remove-stderr", "hide-input", "hide-output"]
using CairoMakie
using GridapMakie
using Printf

using Gridap
using GridapGmsh
using Gridap.Geometry
using GridapPETSc

using Femwell.Maxwell.Electrostatic
using Femwell.Thermal

dir = @__DIR__
run(`python $dir/heater_3d_mesh.py`)

# %% [markdown]
# We can visualize the mesh using GridapMakie:

# %% tags=["remove-stderr", "hide-input", "hide-output"]
function visualize_mesh(model::DiscreteModel)
    volume_tags = ["box","core","heater","via2","via1","metal3","metal2","metal3#e1","metal3#e2"]
    colors = ["grey", "black","orange","green", "green", "brown", "brown", "blue", "blue"]
    alphas = [1.0, 1.0, 0.2, 0.5, 0.5, 0.1, 0.1, 0.5, 0.5]
    fig = Figure(fontsize=20)
    xlims = (-20e-6, 60e-6)
    ylims = (-5e-6, 5e-6)
    zlims = (-5e-6, 5e-6)
    aspect = (xlims[2]-xlims[1])/(ylims[2]-ylims[1])
    tickfunc = values->(x->string(x*1e6)).(values)

    # Mesh
    ∂Ω = BoundaryTriangulation(model, tags=volume_tags)
    ax = fig[1,2] = Axis3(
        fig[1,2], aspect=(aspect, 1.0, 1.0), 
        xlabel="x (μm)", ylabel="\ny (μm)", zlabel="z (μm)\n\n", title="Heater Mesh",
        limits=(xlims, ylims, zlims), titlesize=28,
        xtickformat=tickfunc, ytickformat=tickfunc, ztickformat=tickfunc,
        azimuth=51/40 * pi + pi/20, elevation=pi/8 + pi / 10
    )
    plt = plot!(fig[1,2], ∂Ω, shading=true, title="Mesh", fontsize=20)
    for (i, tag) in enumerate(volume_tags)
        ∂Ω_region = BoundaryTriangulation(model, tags=[tag])
        wireframe!(∂Ω_region, color=(colors[i], alphas[i]), linewidth=0.1, fontsize=20)
    end
    wf = [PolyElement(color=color) for color in colors]
    Legend(fig[1,1], wf, volume_tags, halign=:left, valign=:top, framevisible=false, labelsize=20)

    return fig
end

# %% tags=["remove-stderr"]
model = GmshDiscreteModel("$dir/mesh.msh")
fig = visualize_mesh(model)
display(fig)

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
boundary_temperatures = Dict("metal3#e1___None" => 0.4, "metal3#e2___None" => 0.0)

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

# %% [markdown]
# And we plot the fields using GridapMakie:

# %% tags=["remove-stderr", "hide-input", "hide-output"]
function plot_fields(model, potential, current_density, temperature, heat_flux)

    # Set up figure
    volume_tags = ["box","core","heater","via2","via1","metal3","metal2","metal3#e1","metal3#e2"]
    fontsize = 10*2
    fig = Figure(fontsize=fontsize, resolution=(800*2, 500*2))
    xlims = (-20e-6, 60e-6)
    ylims = (-5e-6, 5e-6)
    zlims = (-5e-6, 5e-6)
    aspect = (xlims[2]-xlims[1])/(ylims[2]-ylims[1])
    tickfunc = values->(x->string(x*1e6)).(values)

    # Set up four quadrants
    g11 = fig[1, 1] = GridLayout()
    g12 = fig[1, 2] = GridLayout()
    g21 = fig[2, 1] = GridLayout()
    g22 = fig[2, 2] = GridLayout()

    # Plot fields
    gs = [g11, g12, g21, g22]
    fs = [potential, current_density, temperature, heat_flux]
    labels = ["Potential", "Current", "Temperature", "Heat Flux"]
    cmaps = [:viridis, :viridis, Reverse(:RdBu), Reverse(:RdBu)]
    cranges = [(0, 0.5), (0, 1.75e10), (0, 120.0), (3.0e4, 4e9)]
    Ω = Triangulation(model, tags=volume_tags)
    for (g, f, label, cmap, crange) ∈ zip(gs, fs, labels, cmaps, cranges)
        Axis3(
            g[1,1], aspect=(aspect, 1.0, 1.0), 
            xlabel="x (μm)", ylabel="\ny (μm)", zlabel="z (μm)\n\n", title=label,
            limits=(xlims, ylims, zlims), titlesize=18*2,
            xtickformat=tickfunc, ytickformat=tickfunc, ztickformat=tickfunc,
            azimuth=51/40 * pi + pi/20, elevation=pi/8 + pi / 10
        )
        plt = plot!(g[1,1], Ω, f, shading=true, colorrange=crange, colormap=cmap, fontsize=fontsize)
        Colorbar(g[1,2], plt, label=label, vertical=true, flipaxis = false)
        colsize!(g, 1, Aspect(1, 1.25))
    end
    return fig
end

# %% tags=["remove-stderr"]
fig = plot_fields(model, potential(p0), current_density(p0), temperature(T0), heat_flux(T0))
display(fig)

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
figure = Figure()
ax = Axis(
    figure[1, 1],
    ylabel = "Average silicon waveguide temperature / K",
    xlabel = "Time / ms",
)

# %% [markdown]
# We can plot the transient 
# %% tags=["remove-stderr", "hide-input", "hide-output"]
function plot_time(uₕₜ::Gridap.ODEs.TransientFETools.TransientFESolution; label="")
    # Set up figure
    volume_tags = ["box","core","heater","via2","via1","metal3","metal2","metal3#e1","metal3#e2"]
    fig = Figure()
    xlims = (-20e-6, 60e-6)
    ylims = (-5e-6, 5e-6)
    zlims = (-5e-6, 5e-6)
    aspect = (xlims[2]-xlims[1])/(ylims[2]-ylims[1])
    tickfunc = values->(x->string(x*1e6)).(values)
    Ω = Triangulation(model, tags=volume_tags)
    fig = Figure()

    # Set up axis
    it = iterate(uₕₜ)
    uh = Any[nothing]
    time = Any[nothing]
    state = Any[nothing]
    (uh[1], time[1]), state[1] = it
    t_uh_time = Observable(time[1])
    timetitle = lift(time -> label * @sprintf("   t = %.2f ms", time * 1e3), t_uh_time)
    Axis3(
        fig[1,1], aspect=(aspect, 1.0, 1.0), 
        xlabel="x (μm)", ylabel="\ny (μm)", zlabel="z (μm)", title=timetitle,
        limits=(xlims, ylims, zlims),
        xtickformat=tickfunc, ytickformat=tickfunc, ztickformat=tickfunc,
        azimuth=51/40 * pi + pi/20, elevation=pi/8 + pi / 10
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
    plt = plot!(fig[1,1], Ω, u, shading=true, colorrange=(0, 110.0), colormap=Reverse(:RdBu))
    Colorbar(fig[1,2], plt, vertical=true, label="Temperature (K)")
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    t_uh_time, fig
end

# %% [markdown]
# Here we set the parameters for the transient simulation
# %% tags=["remove-stderr", "hide-output"]
Δt = 2e-7
t_end = 2e-4
plot_every = 200
timestamps = collect(Δt:Δt*plot_every:t_end)
total_render_time = 5.0
fps = ceil(Int, length(timestamps) / total_render_time)
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
    record(fig, "$dir/animation-$label.gif", timestamps, framerate=fps) do this_t
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
    lines!(ax, t * 1e3, s, label = label)
end

# %% [markdown]
# If we display the animations we can visualize the heatup and cooldown of the device:
# ![heatup animation](animation-heatup.gif)
# ![cooldown animation](animation-cooldown.gif)

# %% [markdown]
# Finally, we can plot an average of the temperature over time: 
# %% tags=["remove-stderr"]
axislegend(position = :rc)
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
