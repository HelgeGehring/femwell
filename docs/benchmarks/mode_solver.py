# ---
# jupyter:
#   jupytext:
#     formats: py:light,md:myst
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# # Benchmark of the mode solver 2

# Reproducing {cite}`Hadley2002`, where the modes of a strip and
# several rib waveguide were calculated and presented with an error value.
# The error for all modes is calculated to be smaller than $\pm 3 \cdot 10^{-6}$,
# thus this should be the maximum derivation for our simulations.
# The results are presented here:

# + tags=["remove-stderr", "hide-input"]
from collections import OrderedDict

import numpy as np
import pandas as pd
import shapely
import shapely.affinity
from juliacall import Main as jl
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

julia_script = [
    "using Gridap",
    "using Gridap.Geometry",
    "using GridapGmsh",
    "using Femwell.Maxwell.Waveguide",
    "using GridapMakie, CairoMakie",
    "import PhysicalConstants.CODATA2018: c_0, μ_0, ε_0",
    "import Unitful: ustrip",
]

for line in julia_script:
    jl.seval(line)

core_width = 3
core_thickness = 1
slab_width = 18
box_thickness = 10
clad_thickness = 3

slab_thicknesses = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
neff_values_paper = [3.412022, 3.412126, 3.412279, 3.412492, 3.412774, 3.413132, 3.413571, 3.414100]
neff_values_femwell_slepc = []
neff_values_femwell_scipy = []
neff_values_julia = []

for slab_thickness in slab_thicknesses:
    if slab_thickness == 0.0:
        core = shapely.Polygon(
            (
                (-core_width / 2, 0),
                (-core_width / 2, core_thickness),
                (core_width / 2, core_thickness),
                (core_width / 2, 0),
            )
        )
    else:
        core = shapely.Polygon(
            (
                (-slab_width / 2, 0),
                (-slab_width / 2, slab_thickness),
                (-core_width / 2, slab_thickness),
                (-core_width / 2, core_thickness),
                (core_width / 2, core_thickness),
                (core_width / 2, slab_thickness),
                (slab_width / 2, slab_thickness),
                (slab_width / 2, 0),
            )
        )

    polygons = OrderedDict(
        core=core,
        box=shapely.box(-slab_width / 2, -box_thickness, slab_width / 2, 0),
        clad=shapely.box(-slab_width / 2, clad_thickness, slab_width / 2, 0),
    )

    resolutions = dict(core={"resolution": 0.1, "distance": 3})

    mesh = from_meshio(
        mesh_from_OrderedDict(
            polygons, resolutions, default_resolution_max=0.3, filename="mesh.msh"
        )
    )
    refractive_indices = {"core": 3.44, "box": 3.40, "clad": 1}
    jl.refractive_indices = refractive_indices

    julia_script = [
        'model = GmshDiscreteModel("mesh.msh")',
        "Ω = Triangulation(model)",
        "labels = get_face_labeling(model)",
        "τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)",
        'refractive_indices = Dict("core" => 3.44^2, "box" => 3.40^2, "clad" => 1^2)',
        "ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in refractive_indices)[tag]",
        "modes = calculate_modes(model, ε ∘ τ,  λ = 1.15, order = 1)",
        "neffs = [real(n_eff(mode)) for mode in modes]",
    ]

    neff_values_julia.append(jl.neffs[0])

    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros()
    for subdomain, n in refractive_indices.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2

    modes = compute_modes(basis0, epsilon, wavelength=1.15, num_modes=1, order=2, solver="slepc")
    neff_values_femwell_slepc.append(np.real(modes[0].n_eff))

    modes = compute_modes(basis0, epsilon, wavelength=1.15, num_modes=1, order=2, solver="scipy")
    neff_values_femwell_scipy.append(np.real(modes[0].n_eff))


pd.DataFrame(
    {
        "slab_thickness": slab_thicknesses,
        "reference value": (f"{n:.6f}" for n in neff_values_paper),
        "calculated value slepc": (f"{n:.6f}" for n in neff_values_femwell_slepc),
        "difference slepc": (
            f"{n1-n2:.6f}" for n1, n2 in zip(neff_values_paper, neff_values_femwell_slepc)
        ),
        "calculated value scipy": (f"{n:.6f}" for n in neff_values_femwell_scipy),
        "difference scipy": (
            f"{n1-n2:.6f}" for n1, n2 in zip(neff_values_paper, neff_values_femwell_scipy)
        ),
        "calculated value julia": (f"{n:.6f}" for n in neff_values_julia),
        "difference julia": (
            f"{n1-n2:.6f}" for n1, n2 in zip(neff_values_paper, neff_values_julia)
        ),
    }
).style.apply(
    lambda differences: [
        "background-color: green" if abs(float(difference)) < 5e-6 else "background-color: red"
        for difference in differences
    ],
    subset=["difference slepc", "difference scipy", "difference julia"],
)
# -

# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
