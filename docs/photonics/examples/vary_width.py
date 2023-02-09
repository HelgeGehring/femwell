# -*- coding: utf-8 -*-
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

# # Variation of the width

# + tags=["remove-stderr"]

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import box
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from tqdm import tqdm

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver import calculate_te_frac, compute_modes

wavelength = 1.55
num_modes = 8
widths = np.linspace(0.5, 3.5, 100)

all_lams = np.zeros((widths.shape[0], num_modes))
all_te_fracs = np.zeros((widths.shape[0], num_modes))
for i, width in enumerate(tqdm(widths)):
    core = box(0, 0, width, 0.5)
    polygons = OrderedDict(
        core=core,
        box=clip_by_rect(core.buffer(1.0, resolution=4), -np.inf, -np.inf, np.inf, 0),
        clad=clip_by_rect(core.buffer(1.0, resolution=4), -np.inf, 0, np.inf, np.inf),
    )

    resolutions = {"core": {"resolution": 0.1, "distance": 1}}

    mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.6))

    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    for subdomain, n in {"core": 1.9963, "box": 1.444, "clad": 1}.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2

    lams, basis, xs = compute_modes(
        basis0, epsilon, wavelength=wavelength, mu_r=1, num_modes=num_modes
    )
    all_lams[i] = np.real(lams)
    all_te_fracs[i, :] = [calculate_te_frac(basis, xs[idx]) for idx in range(num_modes)]

all_lams = np.real(all_lams)
plt.xlabel("Width of waveguide [µm]")
plt.ylabel("Effective refractive index")
plt.ylim(1.444, np.max(all_lams) + 0.1 * (np.max(all_lams) - 1.444))
for lams, te_fracs in zip(all_lams.T, all_te_fracs.T):
    plt.plot(widths, lams)
    plt.scatter(widths, lams, c=te_fracs, cmap="cool")
plt.colorbar().set_label("TE fraction")
plt.show()
