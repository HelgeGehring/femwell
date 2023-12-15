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
# # Effective area and effective refractive index of a Si-NCs waveguide from
# https://opg.optica.org/ol/abstract.cfm?uri=ol-37-12-2295#:~:text=A%20simple%20and%20physically%20meaningful,of%20Si%2DNC%20slot%20waveguides.
# %% tags=["hide-input"]
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import shapely.affinity
from matplotlib.ticker import MultipleLocator

from shapely.ops import clip_by_rect
from femwell.mesh import mesh_from_OrderedDict
from collections import OrderedDict
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from femwell.maxwell.waveguide import compute_modes

# %% [markdown]
# We describe the geometry using shapely.
# In this case it's simple: we use a shapely.box for the waveguide according to the parameter provided in the paper
# For the surrounding we buffer the core and clip it to the part above the waveguide for the air
# Width is a sweep parameter from 50nm to 700nm

# %%
wavelength = 1.55
capital_w = 1.4
capital_h = 0.3
h = 0.5
t = 0.1
n_silicon = 3.48
n_air = 1
n_nc = 1.6
n_silica = 1.45
w_list = [x for x in range(50,700,10)] # x*1e-3 for x in range(50,700,10)
neff_list = []
aeff_list = []
for width in w_list:
    width = width *1e-3
    nc = shapely.geometry.box(-width / 2, capital_h + (h - t) / 2, +width / 2, capital_h + (h - t) / 2 + t)
    silica = shapely.geometry.box(-capital_w / 2, 0, +capital_w / 2, capital_h)
    silicon = shapely.geometry.box(-width / 2, capital_h, +width / 2, capital_h + h)

    polygons = OrderedDict(
        core=nc,
        silicon=silicon,
        silica=silica,
        air=clip_by_rect(nc.buffer(5.0, resolution=5), -np.inf, 0, np.inf, np.inf),
    )

    resolutions = dict(core={"resolution": 0.03, "distance": 0.1})

    mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=10))

    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros()
    for subdomain, n in {"core": n_nc, "silicon": n_silicon, "air": n_air, "silica": n_silica}.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n ** 2
    modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=1, order=2)

    for mode in modes:
        print(f"Effective refractive index: {mode.n_eff:.4f}")
        #print(f"Effective mode area: {mode.calculate_effective_area():.4f}")
        neff_list.append(mode.n_eff)
        #aeff_list.append(mode.calculate_effective_area())
# %% [markdown]
# Plot the result
reference_neff_500nm = pd.read_csv("../reference_data/Rukhlenko/fig_1c_neff/h_500nm.csv",dtype=np.float64)
ref_x,ref_y = np.split(reference_neff_500nm.values, 2, axis=1)

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(9,10))
#ax1.scatter(w_list, aeff_list, c='b',label='stimulated aeff')
#ax1.set_title('aeff at h = 500nm')
ax1.yaxis.set_major_locator(MultipleLocator(0.05))
ax1.yaxis.set_minor_locator(MultipleLocator(0.01))
ax1.set_ylim(0, 0.3)


ax2.scatter(w_list, neff_list, c='r',label='stimulated neff')
ax2.scatter(ref_x,ref_y, c='b',label='reference neff')
ax2.set_title('neff at h = 500nm')
ax2.yaxis.set_major_locator(MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(MultipleLocator(0.2))
ax2.set_ylim(0, 2.8)

plt.legend()
plt.show()