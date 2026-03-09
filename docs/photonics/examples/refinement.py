# ---
# jupyter:
#   jupytext:
#     custom_cell_magics: kql
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: env_3.11
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Mesh refinement

# %% tags=["hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.affinity
from skfem import Basis, ElementTriP0, adaptive_theta
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

# %%
polygons_paper = OrderedDict(
    left=shapely.LineString(((0, 0), (0, 1))),
    right=shapely.LineString(((1, 0), (1, 1))),
    top=shapely.LineString(((0, 1), (1, 1))),
    core=shapely.box(0, 0, 0.5, 0.5),
    clad=shapely.box(0, 0, 1, 1),
)

boundaries = [
    [lambda x: x[0] == np.min(x[0]), lambda x: x[0] == np.max(x[0])],
    [lambda x: x[0] == np.min(x[0]), lambda x: x[0] == np.max(x[0])],
    [lambda x: x[0] == np.min(x[0]), lambda x: x[1] == np.max(x[1])],
    [lambda x: x[0] == np.min(x[0]), lambda x: x[1] == np.max(x[1])],
]

epsilons_paper = [
    {"core": 2.25, "clad": 1},
    {"core": 8, "clad": 1},
    {"core": 1, "clad": 2.25},
    {"core": 1, "clad": 8},
]

neff_values_paper = [1.27627404, 2.65679692, 1.387926425, 2.761465320]
num = 1

mesh = from_meshio(
    mesh_from_OrderedDict(polygons_paper, {}, default_resolution_max=0.08, filename="mesh.msh")
)
mesh.draw().show()
# %%
errors = []
nelements = []

for i in range(20):
    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros()
    for subdomain, e in epsilons_paper[num].items():
        epsilon[basis0.get_dofs(elements=subdomain)] = e

    modes = compute_modes(
        basis0, epsilon, wavelength=1.5, order=1, metallic_boundaries=boundaries[num]
    )
    modes[0].plot(modes[0].E.real, direction="x")
    plt.show()
    error = modes[0].n_eff.real - neff_values_paper[num]
    errors.append(error)
    nelements.append(mesh.nelements)
    elements_to_refine = adaptive_theta(modes[0].eval_error_estimator(), theta=0.5)
    print(
        f"Error {i:2d}: {error: .7f},\tElements: {mesh.nelements:10d}\tSelected: {len(elements_to_refine)}"
    )

    fig, ax = plt.subplots(1, 2)
    for a in ax:
        a.set_aspect("equal")
    mesh.draw(ax=ax[0])
    mesh_restricted = mesh.restrict(elements=elements_to_refine)
    mesh_restricted.draw(color="red", linewidth=2, ax=ax[0])
    mesh = mesh.refined(elements_to_refine)
    mesh.draw(ax=ax[1])
    mesh_restricted.draw(color="red", linewidth=2, ax=ax[1]).show()

fig, axs = plt.subplots(2, 1)
axs[0].plot(errors, label="Error")
axs[0].axhline(y=0, color="r", linestyle="-")
axs[0].set_ylabel("Error")
axs[1].plot(nelements)
axs[1].semilogy()
axs[1].set_xlabel("Iteration")
axs[1].set_ylabel("Number of elements")
plt.show()
# %%
