# %% [markdown]
# # Modes of waveguides with diagonal anisotropy
# Example based on https://ieeexplore.ieee.org/abstract/document/9095209

# %%
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely import box, Polygon, difference
from scipy.constants import epsilon_0, speed_of_light
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0, ElementVector
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict
from femwell.visualization import plot_domains

# %% [markdown]
# We define the geometry and plot the mesh and domains. The waveguide in this example has angled sidewals.

# %%
wg_width = 1.0
wg_thickness = 0.6
slab_thickness = 0.3
sidewall_angle = 60
delta = (wg_thickness - slab_thickness) / 2 / np.tan(np.radians(sidewall_angle))
core = Polygon(
    (
        (wg_width / 2 + delta, wg_thickness - slab_thickness),
        (wg_width / 2 - delta, wg_thickness),
        (-wg_width / 2 + delta, wg_thickness),
        (-wg_width / 2 - delta, wg_thickness - slab_thickness),
    )
)
slab = box(-2, 0, 2, slab_thickness)
env = box(-2, -1, 2, 1.5)

polygons = OrderedDict(
    core=core,
    slab=slab,
    box=clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
)

resolutions = dict(
    core={"resolution": 0.05, "distance": 1.0}, slab={"resolution": 0.1, "distance": 0.5}
)

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=10))
mesh.draw().show()

plot_domains(mesh)
plt.show()

# %% [markdown]
# We define the dielectric constant for Y-propagation, where the extraordinary axis is aligned with "x" in the simulation

# %%
no = 2.2111
ne = 2.1376
ln = np.array([ne, no, no])

epsilons = {"core": ln, "slab": ln, "box": 1.44, "clad": 1.44}

basis0 = Basis(mesh, ElementVector(ElementTriP0(), 3))
epsilon = basis0.zeros()
for subdomain, e in epsilons.items():
    epsilon[basis0.get_dofs(elements=subdomain).all().reshape(-1, 3)] = e**2

(epsilonx, basisx), (epsilony, basisy), (epsilonz, basisz) = basis0.split(epsilon)
basisx.plot(epsilonx, shading="gouraud", colorbar=True).show()
basisy.plot(epsilony, shading="gouraud", colorbar=True).show()
basisz.plot(epsilonz, shading="gouraud", colorbar=True).show()

# %% [markdown]
# Running the simulation, plotting the field and printing the effective indices of the fundamental TE and TM modes

# %%
wavelength = 1.55

modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=2, order=2)
for mode in modes:
    print(f"Effective refractive index: {mode.n_eff:.4f}")
    mode.show("E", part="real", colorbar=True)
    mode.show("E", part="imag", colorbar=True)

# %% [markdown]
# We run again, for Z-propagation

# %%
no = 2.2111
ne = 2.1376
ln = np.array([no, no, ne])

epsilons = {"core": ln, "slab": ln, "box": 1.44, "clad": 1.44}

basis0 = Basis(mesh, ElementVector(ElementTriP0(), 3))
epsilon = basis0.zeros()
for subdomain, e in epsilons.items():
    epsilon[basis0.get_dofs(elements=subdomain).all().reshape(-1, 3)] = e**2

modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=2, order=2)
for mode in modes:
    print(f"Effective refractive index: {mode.n_eff:.4f}")
    mode.show("E", part="real", colorbar=True)
    mode.show("E", part="imag", colorbar=True)
