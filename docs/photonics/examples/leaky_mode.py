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
# # Leaky mode

# %% [markdown]
# Reproducing one example of {cite}`Hu2009`

# %% tags=["remove-stderr", "hide-input"]

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely import LineString, box
from skfem import Basis, ElementDG, ElementTriP1
from skfem.io import from_meshio

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver_2d_periodic import plot_periodic, solve_periodic

# %%

height = 5.76 / 2 + 5
a = 0.010
b = 0.78
c = 0.2
slab = 0.920 + 5
pml = 3

wavelength = 1
k0 = 2 * np.pi / wavelength

left = LineString([(0, y) for y in np.linspace(-height, height, 2)])
right = LineString([(a, y) for y in np.linspace(-height, height, 2)])
top = LineString([(x, height) for x in np.linspace(0, a, 2)])
bottom = LineString([(x, -height) for x in np.linspace(0, a, 2)])

background = box(0, -height, a, height)
structure = box(0, -b / 2, a, b / 2)
structure1 = box(0, height - slab, a, height)
structure2 = box(0, -height + slab, a, -height)

resolutions = {
    "structure": {"resolution": 0.1, "distance": 0.1},
    "hole": {"resolution": 0.1, "distance": 0.1},
}

mesh = from_meshio(
    mesh_from_OrderedDict(
        OrderedDict(
            left=left,
            right=right,
            top=top,
            bottom=bottom,
            structure=structure,
            structure1=structure1,
            structure2=structure2,
            background=background,
        ),
        resolutions=resolutions,
        filename="mesh.msh",
        default_resolution_max=0.05,
        periodic_lines=[("left", "right")],
    )
)

# %%
basis_vec = Basis(mesh, ElementTriP1() * ElementTriP1())
basis_epsilon_r = basis_vec.with_element(ElementDG(ElementTriP1()))

epsilon_r = basis_epsilon_r.zeros(dtype=np.complex64) + 1.45
epsilon_r[basis_epsilon_r.get_dofs(elements="background")] = 1.39
epsilon_r **= 2
basis_epsilon_r.plot(np.real(epsilon_r), aspect=a, colorbar=True).show()

epsilon_r += basis_epsilon_r.project(
    lambda x: (0.5j) * (np.clip(np.abs(x[1]) - height + pml, 0, np.inf) / pml) ** 2,
    dtype=np.complex64,
)
basis_epsilon_r.plot(np.imag(epsilon_r), aspect=a, colorbar=True).show()

ks, basis_phi, phis = solve_periodic(basis_epsilon_r, epsilon_r, k0)

for i, k in enumerate(ks):
    if 0 < np.imag(k) < 0.002:
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        plt.title(f"{k}")
        plot_periodic(k, a, basis_phi, phis[..., i], 100, ax)
        plt.show()

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
