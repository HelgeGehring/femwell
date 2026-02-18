module Thermal

using Gridap
using Gridap.Algebra
using GridapODEs.ODETools
using GridapODEs.TransientFETools

export calculate_temperature, temperature, calculate_temperature_transient, heat_flux

struct Temperature
    temperature::FEFunction
    thermal_conductivity::CellField
end

function calculate_temperature(
    thermal_conductivity::CellField,
    power_density::CellField,
    temperatures = Dict{String,Float64};
    order::Int = 1,
    solver::LinearSolver = BackslashSolver(),
)
    tags = collect(keys(temperatures))
    tags_temperatures = [temperatures[tag] for tag in tags]

    model = get_active_model(get_triangulation(thermal_conductivity))
    V = TestFESpace(
        model,
        ReferenceFE(lagrangian, Float64, order);
        conformity = :H1,
        dirichlet_tags = tags,
    )
    U = TrialFESpace(V, tags_temperatures)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, order)

    a(u, v) = ∫(∇(v) ⋅ thermal_conductivity ⋅ ∇(u))dΩ
    b(v) = ∫(v * power_density)dΩ

    op = AffineFEOperator(a, b, U, V)
    temperature = solve(solver, op)

    return Temperature(temperature, thermal_conductivity)
end

temperature(temperature::Temperature) = temperature.temperature
heat_flux(temperature::Temperature) =
    temperature.thermal_conductivity * ∇(temperature.temperature)

function calculate_temperature_transient(
    thermal_conductivity::CellField,
    thermal_diffusitivity::CellField,
    power_density::Union{Function,CellField,Float64},
    temperatures::Dict{String,Float64},
    T0::CellField,
    Δt::Float64,
    t_end::Float64;
    order::Int = 1,
    solver::LinearSolver = LUSolver(),
)
    temperatures_c = collect(temperatures)
    tags = [u for (u, v) in temperatures_c]

    function def_g(value)
        g(x, t::Real) = value
        g(t::Real) = x -> g(x, t)
        return g
    end

    tags_temperatures = [def_g(v) for (u, v) in temperatures_c]
    if isempty(tags_temperatures)
        tags_temperatures = x -> 0
    end

    model = get_active_model(get_triangulation(thermal_diffusitivity))
    V = TestFESpace(
        model,
        ReferenceFE(lagrangian, Float64, order);
        conformity = :H1,
        dirichlet_tags = tags,
    )
    U = TransientTrialFESpace(V, tags_temperatures)


    Ω = Triangulation(model)
    dΩ = Measure(Ω, order)
    m₀(u, v) = ∫(thermal_conductivity / thermal_diffusitivity * u * v)dΩ
    a₀(u, v) = ∫(thermal_conductivity * (∇(u) ⋅ ∇(v)))dΩ

    op_C = if typeof(power_density) <: Function
        println("function!")
        TransientConstantMatrixFEOperator(m₀, a₀, (t, v) -> ∫(power_density(t) * v)dΩ, U, V)
    else
        TransientConstantFEOperator(m₀, a₀, v -> ∫(power_density * v)dΩ, U, V)
    end

    θ = 0.5
    ode_solver = ThetaMethod(solver, Δt, θ)

    u₀ = interpolate_everywhere(T0, U(0.0))
    uₕₜ = solve(ode_solver, op_C, u₀, 0, t_end)

    return uₕₜ
end

end
