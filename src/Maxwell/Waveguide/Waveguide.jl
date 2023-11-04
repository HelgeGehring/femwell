module Waveguide

using Gridap
using Gridap.Geometry
using Gridap.TensorValues
using Arpack
using GridapMakie, CairoMakie
import PhysicalConstants.CODATA2018: c_0, μ_0, ε_0
import Unitful: ustrip

export Mode
export calculate_modes
export E, H, power, overlap, te_fraction, tm_fraction, coupling_coefficient, perturbed_neff
export frequency, n_eff, ω, λ
export plot_mode, plot_field, write_mode_to_vtk

struct Mode
    k0::ComplexF64
    k::ComplexF64
    E::FEFunction
    order::Int
    ε::CellField
end

Base.show(io::IO, mode::Mode) = print(io, "Mode(k=$(mode.k), n_eff=$(n_eff(mode)))")
Gridap.Measure(mode::Mode) = Measure(get_triangulation(mode.E), 2 * mode.order)
n_eff(mode::Mode) = mode.k / mode.k0
ω(mode::Mode) = mode.k0 * ustrip(c_0)
frequency(mode::Mode) = 2π * ω(mode)
λ(mode::Mode) = 2π / mode.k
E(mode::Mode) =
    mode.E[1] ⋅ TensorValue([1.0 0.0 0.0; 0.0 1.0 0.0]) +
    mode.E[2] * VectorValue(0.0, 0.0, 1.0)
H(mode::Mode) =
    -1im / ustrip(μ_0) / ω(mode) * (
        (1im * mode.k * mode.E[1] - ∇(mode.E[2])) ⋅
        TensorValue([0.0 1.0 0.0; -1.0 0.0 0.0]) +
        ∇(mode.E[1]) ⊙ TensorValue(0.0, -1.0, 1.0, 0.0) * VectorValue(0.0, 0.0, 1.0)
    )
overlap(mode1::Mode, mode2::Mode) =
    0.5 * sum(
        ∫(
            (cross(conj(E(mode1)), H(mode2)) + cross(E(mode2), conj(H(mode1)))) ⋅
            VectorValue(0.0, 0.0, 1.0),
        )Measure(mode1),
    )
coupling_coefficient(mode1::Mode, mode2::Mode, delta_epsilon) =
    sum(∫(delta_epsilon * conj(E(mode1)) ⋅ E(mode2))Measure(mode1))
perturbed_neff(mode::Mode, delta_epsilon::CellField) =
    n_eff(mode) + coupling_coefficient(mode, mode, delta_epsilon) * ustrip(ε_0 * c_0) * 0.5
intensity(mode::Mode) = real(
    (cross(conj(E(mode)), H(mode)) + cross(E(mode), conj(H(mode)))) ⋅ VectorValue(0, 0, 1),
)
power(mode::Mode) = overlap(mode, mode)
propagation_loss(mode::Mode, distance::Real) =
    -20 / log(10) * mode.k0 * np.imag(n_eff(mode)) * distance

te_fraction(mode::Mode) = abs(
    ∑(∫(mode.E[1] ⋅ TensorValue(1, 0, 0, 0) ⋅ mode.E[1])Measure(mode)) /
    ∑(∫(mode.E[1] ⋅ mode.E[1])Measure(mode)),
)
tm_fraction(mode::Mode) = abs(
    ∑(∫(mode.E[1] ⋅ TensorValue(0, 0, 0, 1) ⋅ mode.E[1])Measure(mode)) /
    ∑(∫(mode.E[1] ⋅ mode.E[1])Measure(mode)),
)

function calculate_modes(
    model::DiscreteModel,
    ε::CellField;
    k0::Union{Number,Nothing} = nothing,
    λ::Union{Number,Nothing} = nothing,
    order::Int = 1,
    num::Int = 1,
    radius::Real = Inf,
    k0_guess::Union{Number,Nothing} = nothing,
    metallic_boundaries = String[],
)
    if count(isnothing, [k0, λ]) != 1
        throw(ArgumentError("Exactly one of k0,λ must be defined"))
    end
    k0 = !isnothing(k0) ? k0 : 2π / λ

    V1 = TestFESpace(
        model,
        ReferenceFE(nedelec, order),
        vector_type = Vector{ComplexF64},
        dirichlet_tags = metallic_boundaries,
    )
    U1 = TrialFESpace(V1)
    V2 = TestFESpace(
        model,
        ReferenceFE(lagrangian, Float64, order + 1),
        vector_type = Vector{ComplexF64},
        dirichlet_tags = metallic_boundaries,
    )
    U2 = TrialFESpace(V2)
    V = MultiFieldFESpace([V1, V2])
    U = MultiFieldFESpace([U1, U2])
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2 * order)

    epsilons = ε(get_cell_points(Measure(Ω, 1)))
    k0_guess =
        isnothing(k0_guess) ? k0^2 * maximum(maximum.(maximum.(real.(epsilons)))) * 1.1 :
        k0_guess

    lhs, rhs = if radius == Inf
        μ_r = 1
        twothree = TensorValue([1 0 0; 0 1 0])
        lhs_straight((u1, u2), (v1, v2)) =
            ∫(
                1 / μ_r * (curl(u1) ⊙ curl(v1) / k0^2 + ∇(u2) ⊙ v1) + (
                    -u1 ⋅ twothree ⋅ ε ⋅ transpose(twothree) ⋅ v1 +
                    u1 ⋅ twothree ⋅ ε ⋅ transpose(twothree) ⋅ ∇(v2) -
                    u2 * VectorValue(0, 0, 1) ⋅ ε ⋅ VectorValue(0, 0, 1) * v2 * k0^2
                ),
            )dΩ
        rhs_straight((u1, u2), (v1, v2)) = ∫(-1 / μ_r * u1 ⊙ v1 / k0^2)dΩ

        lhs_straight, rhs_straight
    else
        radius_function(x) = x[1]
        r = interpolate_everywhere(radius_function, V2)
        lhs_radial((u1, u2), (v1, v2)) =
            ∫(
                -(r / radius * curl(u1) ⋅ curl(v1)) +
                (radius / r * gradient(u2) ⋅ v1) +
                (k0^2 * ε * r / radius * u1 ⋅ v1) -
                (radius / r * gradient(u2) ⋅ gradient(v2)) +
                (k0^2 * ε * radius / r * u2 * v2),
            )dΩ

        rhs_radial((u1, u2), (v1, v2)) =
            ∫((radius / r * u1 ⋅ v1) - (radius / r * u1 ⋅ gradient(v2)))dΩ

        lhs_radial, rhs_radial
    end


    assem = Gridap.FESpaces.SparseMatrixAssembler(U, V)
    A = assemble_matrix(lhs, assem, U, V)
    B = assemble_matrix(rhs, assem, U, V)
    if all(imag(A.nzval) .== 0) && all(imag(B.nzval) .== 0)
        A, B = real(A), real(B)
    end
    vals, vecs = eigs(A, B, sigma = k0_guess, nev = num)

    if radius == Inf
        vecs[num_free_dofs(V1)+1:end, :] ./= 1im * sqrt.(vals)' / k0^2
    else
        vecs[num_free_dofs(V1)+1:end, :] ./= 1im * sqrt.(vals)'
    end


    return [
        Mode(
            k0,
            sqrt(k2),
            FEFunction(U, E / sqrt(power(Mode(k0, sqrt(k2), FEFunction(U, E), order, ε)))),
            order,
            ε,
        ) for (k2, E) in zip(vals, eachcol(vecs))
    ]
end

function plot_field(field)
    fig, _, plt = plot(get_triangulation(field), real(field))
    Colorbar(fig[1, 2], plt)
    display(fig)
end

function plot_mode(
    mode::Mode;
    vertical = false,
    vectors = false,
    same_range = false,
    absolute = false,
)
    Ω = get_triangulation(mode.E)
    model = get_active_model(Ω)
    labels = get_face_labeling(model)
    boundary_tags =
        setdiff(unique(get_face_tag(labels, 1)), unique(get_face_tag(labels, 2)))
    ∂Ω = BoundaryTriangulation(model, tags = boundary_tags)

    fig = Figure()
    if !vectors
        fields = [
            ("E_x", VectorValue(1, 0, 0)),
            ("E_y", VectorValue(0, 1, 0)),
            ("E_z", VectorValue(0, 0, 1im)),
        ]

        colorranges = Float64[]
        for (i, (title, vector)) in enumerate(fields)
            efield = real((E(mode) ⋅ vector)(get_cell_points(Measure(Ω, 1))))
            colorrange = max(abs(maximum(maximum.(efield))), abs(minimum(minimum.(efield))))
            push!(colorranges, colorrange)
        end

        for (i, (title, vector)) in enumerate(fields)
            colorrange = same_range ? maximum(colorranges) : colorranges[i]

            x = !vertical + i * vertical
            y = vertical + i * !vertical

            ax = Axis(fig[x, y], title = title)
            plt = plot!(
                ax,
                Ω,
                (absolute ? abs : real)((E(mode) ⋅ vector)),
                #colorrange = (-colorrange, colorrange),
            )
            wireframe!(fig[x, y], ∂Ω, color = :black)
            Colorbar(fig[x+!vertical, y+vertical], plt, vertical = vertical)
        end
    else
        cellpoints = get_cell_points(Measure(Ω, 1))
        xmin = minimum(map(x -> x.data[1], getindex.(cellpoints.cell_phys_point, 1)))
        xmax = maximum(map(x -> x.data[1], getindex.(cellpoints.cell_phys_point, 1)))
        ymin = minimum(map(x -> x.data[2], getindex.(cellpoints.cell_phys_point, 1)))
        ymax = maximum(map(x -> x.data[2], getindex.(cellpoints.cell_phys_point, 1)))

        streamplot(
            (x, y) ->
                Point2(collect(real(∇(mode.E[1] ⋅ mode.E[1]))(Gridap.Point(x, y)).data)),
            xmin .. xmax,
            ymin .. ymax,
        )

        efield = real((E(mode) ⋅ VectorValue(0, 0, 1im))(get_cell_points(Measure(Ω, 1))))
        colorrange = max(abs(maximum(maximum.(efield))), abs(minimum(minimum.(efield))))

        ax = Axis(fig[1, 2], title = "E_z")
        plt = plot!(
            ax,
            Ω,
            real((E(mode) ⋅ VectorValue(0, 0, 1im))),
            colorrange = (-colorrange, colorrange),
        )
        wireframe!(fig[1, 2], ∂Ω, color = :black)
        Colorbar(fig[2, 2], plt, vertical = false)
    end
    display(fig)
end


Base.real(x::VectorValue{D,ComplexF64}) where {D} = VectorValue(real.(x.data))
Base.imag(x::VectorValue{D,ComplexF64}) where {D} = VectorValue(imag.(x.data))
Base.real(x::TensorValue) = TensorValue(real.(x.data))
Base.imag(x::TensorValue) = TensorValue(imag.(x.data))

function write_mode_to_vtk(filename::String, mode::Mode)
    writevtk(
        get_triangulation(mode.E),
        filename,
        cellfields = [
            "E_real" => real(E(mode)),
            "E_imag" => imag(E(mode)),
            "H_real" => real(H(mode)),
            "H_imag" => imag(H(mode)),
            "Intensity" => intensity(mode),
            "ε_real" => real(mode.ε),
            "ε_imag" => imag(mode.ε),
        ],
    )
end

# https://www.nature.com/articles/ncomms8027#Sec11
A1(mode::Mode) = ∑(
    ∫(
        ω(mode) * E(mode) ⋅ TensorValue([1 0 0; 0 1 0; 0 0 0]) ⋅ (ustrip(ε_0) * mode.ε) ⋅ E(mode) -
        1 / ω(mode) * curl(mode.E[1]) ⋅ 1 / ustrip(μ_0) ⋅ curl(mode.E[1]),
    )Measure(mode),
)
A2(mode::Mode) = ∑(
    ∫(
        ω(mode) * ustrip(μ_0) * H(mode) ⋅ TensorValue([1 0 0; 0 1 0; 0 0 0]) ⋅ H(mode) +
        1 / ω(mode) * mode.E[2] ⋅ (ustrip(ε_0) * mode.ε) ⋅ mode.E[2] * ω(mode)^2,
    )Measure(mode),
)
B(mode::Mode) = 2∑(∫(E(mode) ⋅ TensorValue([0 1 0; 0 -1 0; 0 0 0]) ⋅ H(mode))Measure(mode))

end
