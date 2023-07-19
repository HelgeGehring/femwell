module Thermal

using Gridap

export calculate_temperature, temperature

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
    b(v) = ∫(v * conductivity * power_density)dΩ

    op = AffineFEOperator(a, b, U, V)
    temperature = solve(op)

    return Temperature(temperature)
end

temperature(temperature::Temperature) = temperature.temperature

end
