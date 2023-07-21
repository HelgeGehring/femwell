module Thermal

using Gridap

export calculate_temperature, temperature, calculate_temperature_transient

struct Temperature
    temperature::FEFunction
end

function calculate_temperature(
    conductivity::CellField,
    power_density::CellField,
    temperatures = Dict{String,Float64},
    order::Int = 1,
)
    tags = collect(keys(temperatures))
    tags_temperatures = [temperatures[tag] for tag in tags]

    model = get_active_model(get_triangulation(conductivity))
    V = TestFESpace(
        model,
        ReferenceFE(lagrangian, Float64, order);
        conformity = :H1,
        dirichlet_tags = tags,
    )
    U = TrialFESpace(V, tags_temperatures)
    Ω = Triangulation(model)
    dΩ = Measure(Ω, order)

    a(u, v) = ∫(conductivity * ∇(v) ⋅ ∇(u))dΩ
    b(v) = ∫(v * power_density)dΩ

    op = AffineFEOperator(a, b, U, V)
    temperature = solve(op)

    return Temperature(temperature)
end

temperature(temperature::Temperature) = temperature.temperature

function calculate_temperature_transient(
    conductivity::CellField,
    diffusitivity::CellField,
    power_density::Union{CellField,Float64},
    temperatures::Dict{String,Float64},
    T0::CellField,
    Δt::Float64,
    t_end::Float64,
    order::Int = 1,
)
    temperatures_c = collect(temperatures)
    tags = [u for (u, v) in temperatures_c]

    function def_g(value)
        g(x, t::Real) = value
        g(t::Real) = x -> g(x, t)
        return g
    end

    tags_temperatures = [def_g(v) for (u, v) in temperatures_c]

    model = get_active_model(get_triangulation(diffusitivity))
    V = TestFESpace(
        model,
        ReferenceFE(lagrangian, Float64, order);
        conformity = :H1,
        dirichlet_tags = tags,
    )
    U = TransientTrialFESpace(V, tags_temperatures)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, order)
    m₀(u, v) = ∫(conductivity / diffusitivity * u * v)dΩ
    b₀(v) = ∫(power_density * v)dΩ
    a₀(u, v) = ∫(conductivity * (∇(u) ⋅ ∇(v)))dΩ
    op_C = TransientConstantFEOperator(m₀, a₀, b₀, U, V)

    linear_solver = LUSolver()
    θ = 0.5
    ode_solver = ThetaMethod(linear_solver, Δt, θ)

    u₀ = interpolate_everywhere(T0, U(0.0))
    uₕₜ = solve(ode_solver, op_C, u₀, 0, t_end)

    return uₕₜ
end

end
