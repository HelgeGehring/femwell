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
# # Effective area and effective refractive index of a Si-NCs waveguide from {cite}`Rukhlenko2012`.
# %% tags=["hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely.affinity
from matplotlib.ticker import MultipleLocator
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

# %% [markdown]
# We describe the geometry using shapely.
# In this case it's simple: we use a shapely.box for the waveguide according to the parameter provided in the paper
# For the surrounding we buffer the core and clip it to the part above the waveguide for the air
# Width is a sweep parameter from 50nm to 700nm

# %%  tags=["hide-output"]
wavelength = 1.55
capital_w = 1.4
capital_h = 0.3
h_list = [0.5, 0.7]
t = 0.1
n_silicon = 3.48
n_air = 1
n_nc = 1.6
n_silica = 1.45
w_list = [x for x in range(250, 700, 100)]
neff_dict = dict()
aeff_dict = dict()
tm_dict = dict()
for h in h_list:
    neff_list = []
    aeff_list = []
    tm_list = []
    for width in w_list:
        width = width * 1e-3
        nc = shapely.geometry.box(
            -width / 2, capital_h + (h - t) / 2, +width / 2, capital_h + (h - t) / 2 + t
        )
        silica = shapely.geometry.box(-capital_w / 2, 0, +capital_w / 2, capital_h)
        silicon = shapely.geometry.box(-width / 2, capital_h, +width / 2, capital_h + h)

        polygons = OrderedDict(
            core=nc,
            silicon=silicon,
            silica=silica,
            air=nc.buffer(10.0, resolution=4),
        )

        resolutions = dict(
            core={"resolution": 0.02, "distance": 0.3},
            silicon={"resolution": 0.02, "distance": 0.1},
            silica={"resolution": 0.04, "distance": 0.5},
        )

        mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=2))

        basis0 = Basis(mesh, ElementTriP0())
        epsilon = basis0.zeros()
        for subdomain, n in {
            "core": n_nc,
            "silicon": n_silicon,
            "air": n_air,
            "silica": n_silica,
        }.items():
            epsilon[basis0.get_dofs(elements=subdomain)] = n**2
        modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)

        for mode in modes:
            if mode.tm_fraction > 0.5:
                # mode.show(np.real(mode.E))
                print(f"Effective refractive index: {mode.n_eff:.4f}")
                print(f"Effective mode area: {mode.calculate_effective_area(field='y'):.4f}")
                print(f"Mode transversality: {mode.transversality}")
                neff_list.append(np.real(mode.n_eff))
                aeff_list.append(mode.calculate_effective_area())
                tm_list.append(mode.transversality)
                break
        else:
            print(f"no TM mode found for {width}")
    neff_dict[str(h)] = neff_list
    aeff_dict[str(h)] = aeff_list
    tm_dict[str(h)] = tm_list
# %% [markdown]
# Plot the result

# %%
reference_neff_500nm = pd.read_csv(
    "../reference_data/Rukhlenko/fig_1c_neff/h_500nm.csv", dtype=np.float64
)
reference_aeff_500nm = pd.read_csv(
    "../reference_data/Rukhlenko/fig_1b_aeff/0.5_Eq2.csv", dtype=np.float64
)
reference_tm_500nm = pd.read_csv(
    "../reference_data/Rukhlenko/fig_1c_neff/tm_h_500nm.csv", dtype=np.float64
)

reference_neff_700nm = pd.read_csv(
    "../reference_data/Rukhlenko/fig_1c_neff/h_700nm.csv", dtype=np.float64
)
reference_aeff_700nm = pd.read_csv(
    "../reference_data/Rukhlenko/fig_1b_aeff/0.7 Eq2.csv", dtype=np.float64
)
reference_tm_700nm = pd.read_csv(
    "../reference_data/Rukhlenko/fig_1c_neff/tm_h_700nm.csv", dtype=np.float64
)

ref_neff_500nm_x, ref_neff_500nm_y = np.split(reference_neff_500nm.values, 2, axis=1)
ref_aeff_500nm_x, ref_aeff_500nm_y = np.split(reference_aeff_500nm.values, 2, axis=1)
ref_tm_500nm_x, ref_tm_500nm_y = np.split(reference_tm_500nm.values, 2, axis=1)

ref_neff_700nm_x, ref_neff_700nm_y = np.split(reference_neff_700nm.values, 2, axis=1)
ref_aeff_700nm_x, ref_aeff_700nm_y = np.split(reference_aeff_700nm.values, 2, axis=1)
ref_tm_700nm_x, ref_tm_700nm_y = np.split(reference_tm_700nm.values, 2, axis=1)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(9, 20))

ax1.plot(w_list, aeff_dict["0.5"], "-o", c="r", label="stimulated aeff")
ax1.scatter(ref_aeff_500nm_x, ref_aeff_500nm_y, c="b", label="reference aeff")
ax1.set_title("aeff at h = 500nm")
ax1.yaxis.set_major_locator(MultipleLocator(0.05))
ax1.yaxis.set_minor_locator(MultipleLocator(0.01))
ax1.set_ylim(0, 0.3)
ax1.legend()

ax2b = ax2.twinx()
ax2b.set_ylabel("mode transversality")
ax2b.scatter(ref_tm_500nm_x, ref_tm_500nm_y, marker="v", c="b", label="reference transversality")
ax2b.plot(w_list, tm_dict["0.5"], "-v", c="r", label="stimulated transversality")
ax2b.set_ylim(0.775, 1)
ax2b.legend()

ax2.plot(w_list, neff_dict["0.5"], "-o", c="r", label="stimulated neff")
ax2.scatter(ref_neff_500nm_x, ref_neff_500nm_y, c="b", label="reference neff")
ax2.set_title("neff and mode transversality at h = 500nm")
ax2.set_ylabel("neff")
ax2.yaxis.set_major_locator(MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(MultipleLocator(0.2))
ax2.set_ylim(0, 2.8)
ax2.legend()

ax3.plot(w_list, aeff_dict["0.7"], "-o", c="r", label="stimulated aeff")
ax3.scatter(ref_aeff_700nm_x, ref_aeff_700nm_y, c="b", label="reference aeff")
ax3.set_title("aeff at h = 700nm")
ax3.yaxis.set_major_locator(MultipleLocator(0.05))
ax3.yaxis.set_minor_locator(MultipleLocator(0.01))
ax3.set_ylim(0, 0.3)
ax3.legend()

ax4b = ax4.twinx()
ax4b.set_ylabel("mode transversality")
ax4b.scatter(ref_tm_700nm_x, ref_tm_700nm_y, marker="v", c="b", label="reference transversality")
ax4b.plot(w_list, tm_dict["0.7"], "-v", c="r", label="stimulated transversality")
ax4b.set_ylim(0.775, 1)
ax4b.legend()


ax4.plot(w_list, neff_dict["0.7"], "-o", c="r", label="stimulated aeff")
ax4.scatter(ref_neff_700nm_x, ref_neff_700nm_y, c="b", label="reference aeff")
ax4.set_title("neff and mode transversality at h = 700nm")
ax2.yaxis.set_major_locator(MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(MultipleLocator(0.2))
ax4.set_ylim(0, 2.8)
ax4.legend()

plt.show()
