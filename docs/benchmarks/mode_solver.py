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

# # Benchmark of the mode solver

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
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver import compute_modes

core_width = 3
core_thickness = 1
slab_width = 18
box_thickness = 10
clad_thickness = 3

slab_thicknesses = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
neff_values_paper = [3.412022, 3.412126, 3.412279, 3.412492, 3.412774, 3.413132, 3.413571, 3.414100]
neff_values_femwell = []

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

    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros()
    for subdomain, n in {"core": 3.44, "box": 3.40, "clad": 1}.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2

    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.15, num_modes=1, order=2)

    neff_values_femwell.append(np.real(lams[0]))

pd.DataFrame(
    {
        "slab_thickness": slab_thicknesses,
        "reference value": (f"{n:.6f}" for n in neff_values_paper),
        "calculated value": (f"{n:.6f}" for n in neff_values_femwell),
        "difference": (f"{n1-n2:.6f}" for n1, n2 in zip(neff_values_paper, neff_values_femwell)),
    }
).style.apply(
    lambda differences: [
        "background-color: green" if 3 < 1e-6 else "background-color: red" for i in differences
    ],
    subset=["difference"],
)
# -

# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
