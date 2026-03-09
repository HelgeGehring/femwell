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
# # Coax cable

# In this example we calculate the propagation constant and the characteristic impedance of a RG-6 coaxial cable.
# The cable has an 18 AWG (1.024 mm) center conductor and is designed to have a characteristic impedance of $75\Omega$.

# %% tags = ["hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import shapely
import shapely.ops
from skfem import Basis, ElementTriP0, Functional
from skfem.helpers import inner
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict


# %%
def mesh_coax(filename, radius_inner, radius_outer):
    core = shapely.Point(0, 0).buffer(radius_inner)
    isolator2 = shapely.Point(0, 0).buffer((radius_outer + radius_inner) / 4, resolution=32)
    isolator = shapely.Point(0, 0).buffer(radius_outer)

    polygons = OrderedDict(
        surface=shapely.LineString(isolator.exterior),
        surface_core=shapely.LinearRing(core.exterior),
        isolator2=isolator2 - core,
        isolator=isolator - core,
    )

    resolutions = dict(
        isolator2={"resolution": 0.05, "distance": 1},
        isolator={"resolution": 0.1, "distance": 1},
    )

    return mesh_from_OrderedDict(
        polygons, resolutions, filename=filename, default_resolution_max=1, mesh_scaling_factor=1e-3
    )


# %%
frequency = 10e9

mesh = from_meshio(mesh_coax(radius_inner=0.512, radius_outer=2.23039, filename="mesh.msh"))

basis0 = Basis(mesh, ElementTriP0(), intorder=4)
epsilon = basis0.zeros() + 1.29
basis0.plot(epsilon, colorbar=True).show()

conductors = ["isolator2___isolator"]
modes = compute_modes(
    basis0,
    epsilon,
    wavelength=scipy.constants.speed_of_light / frequency,
    mu_r=1,
    num_modes=len(conductors),
    metallic_boundaries=("surface", "surface_core"),
)
print("Propagation constant", 1 / modes.n_effs)

modes[0].show("E", part="real", plot_vectors=True, colorbar=True)

# %%


@Functional(dtype=np.complex64)
def current_form(w):
    return inner(np.array([w.n[1], -w.n[0]]), w.H)


currents = np.zeros((len(conductors), len(modes)))

for mode_i, mode in enumerate(modes):
    modes[0].show("H", part="real", plot_vectors=True, colorbar=True)

    (ht, ht_basis), (hz, hz_basis) = mode.basis.split(mode.H)
    for conductors_i, conductor in enumerate(conductors):
        facet_basis = ht_basis.boundary(facets=mesh.boundaries[conductor])
        current = abs(current_form.assemble(facet_basis, H=facet_basis.interpolate(ht)))
        currents[conductors_i, mode_i] = current

characteristic_impedances = np.linalg.inv(currents).T @ np.linalg.inv(currents)
print("Characteristic impedances", characteristic_impedances)

# %%
