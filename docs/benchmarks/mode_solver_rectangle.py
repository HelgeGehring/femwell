# ---
# jupyter:
#   jupytext:
#     formats: py:percent,md:myst
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# %% [markdown]
# # Benchmark of the mode solver 1

# %% [markdown]
# Reproducing {cite}`Hadley2002`, where the modes of a analytically solvable geometry are calculated.
# The error for all modes is calculated to be smaller than $\pm 1 \cdot 10^{-8}$.
# We'll show that we get pretty close, but will stop at a resonable resolution to keep the runtime sensible.
# Getting even higher accurancy will be left open for adaptive refinement.
# The results are presented here:

# %% tags=["remove-stderr", "hide-input"]
from collections import OrderedDict

import numpy as np
import pandas as pd
import shapely
import shapely.affinity
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

epsilons_paper = [
    {"core": 2.25, "clad": 1},
    {"core": 8, "clad": 1},
    {"core": 1, "clad": 2.25},
    {"core": 1, "clad": 8},
]
boundaries = [["left", "right"], ["left", "right"], ["left", "top"], ["left", "top"]]
neff_values_paper = [1.27627404, 2.65679692, 1.387926425, 2.761465320]
neff_values_femwell_slepc = []
neff_values_femwell_scipy = []

polygons = OrderedDict(
    left=shapely.LineString(((0, 0), (0, 1))),
    right=shapely.LineString(((1, 0), (1, 1))),
    top=shapely.LineString(((0, 1), (1, 1))),
    core=shapely.box(0, 0, 0.5, 0.5),
    clad=shapely.box(0, 0, 1, 1),
)

mesh = from_meshio(
    mesh_from_OrderedDict(polygons, {}, default_resolution_max=0.01, filename="mesh.msh")
)

basis0 = Basis(mesh, ElementTriP0())

for epsilons, boundaries in zip(epsilons_paper, boundaries):
    epsilon = basis0.zeros()
    for subdomain, e in epsilons.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = e

    modes = compute_modes(
        basis0,
        epsilon,
        wavelength=1.5,
        num_modes=1,
        order=2,
        metallic_boundaries=boundaries,
        solver="slepc",
    )
    neff_values_femwell_slepc.append(np.real(modes[0].n_eff))

    modes = compute_modes(
        basis0,
        epsilon,
        wavelength=1.5,
        num_modes=1,
        order=2,
        metallic_boundaries=boundaries,
        solver="scipy",
    )
    neff_values_femwell_scipy.append(np.real(modes[0].n_eff))

pd.DataFrame(
    {
        "Epsilons": [
            f"{epsilons['core']:.2f} / {epsilons['clad']:.2f}" for epsilons in epsilons_paper
        ],
        "Reference value": (f"{n:.8f}" for n in neff_values_paper),
        "Calculated value slepc": (f"{n:.8f}" for n in neff_values_femwell_slepc),
        "Difference slepc": (
            f"{n1-n2:.8f}" for n1, n2 in zip(neff_values_paper, neff_values_femwell_slepc)
        ),
        "Calculated value scipy": (f"{n:.8f}" for n in neff_values_femwell_scipy),
        "Difference scipy": (
            f"{n1-n2:.8f}" for n1, n2 in zip(neff_values_paper, neff_values_femwell_scipy)
        ),
    }
).style.apply(
    lambda differences: [
        "background-color: green" if abs(float(difference)) < 5e-6 else "background-color: red"
        for difference in differences
    ],
    subset=["difference slepc", "difference scipy"],
)

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
