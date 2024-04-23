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
# # Crosstalk reduction
# This example reproduces Figure 3 from {cite}`Yang2020`` in which the crosstalk between adjacent waveguides is suppressed.

# %% tags=["hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon, box
from skfem import Basis, ElementTriP0
from skfem.io import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

# %%
w_sim = 4
h_clad = 1
h_box = 1
w_core_1 = 0.5
w_core_2 = 0.5
gap = 0.5
h_core = 0.22
w_core_c = 0.2
num_c = 2

wavelength = 1.55

references = ("3a.csv", "3b.csv")

for num_c, reference in zip([1, 2], references):
    reference_data = np.loadtxt(reference, unpack=True, delimiter=",")

    w_core_cs = np.linspace(np.min(reference_data[0]), np.max(reference_data[0]), 15) * 1e-3
    coupling_lengths = []

    for w_core_c in w_core_cs:

        def core_c_pos(i):
            return (
                -gap / 2 + (gap - w_core_c * num_c) * (i + 1) / (num_c + 1) + (i + 0.5) * w_core_c
            )

        polygons = OrderedDict(
            core_1=Polygon(
                [
                    (-w_core_1 - gap / 2, 0),
                    (-w_core_1 - gap / 2, h_core),
                    (-gap / 2, h_core),
                    (-gap / 2, 0),
                ]
            ),
            core_2=Polygon(
                [
                    (w_core_2 + gap / 2, 0),
                    (w_core_2 + gap / 2, h_core),
                    (gap / 2, h_core),
                    (gap / 2, 0),
                ]
            ),
            **{
                f"core_c_{i}": Polygon(
                    [
                        (core_c_pos(i) - w_core_c / 2, 0),
                        (core_c_pos(i) - w_core_c / 2, h_core),
                        (core_c_pos(i) + w_core_c / 2, h_core),
                        (core_c_pos(i) + w_core_c / 2, 0),
                    ]
                )
                for i in range(num_c)
            },
            clad=Polygon(
                [(-w_sim / 2, 0), (-w_sim / 2, h_clad), (w_sim / 2, h_clad), (w_sim / 2, 0)]
            ),
            box=Polygon(
                [(-w_sim / 2, 0), (-w_sim / 2, -h_box), (w_sim / 2, -h_box), (w_sim / 2, 0)]
            ),
        )

        resolutions = dict(
            core_1={"resolution": 0.02, "distance": 0.3, "SizeMax": 0.1},
            core_2={"resolution": 0.02, "distance": 0.3, "SizeMax": 0.1},
            **{
                f"core_c_{i}": {"resolution": 0.01, "distance": 0.3, "SizeMax": 0.1}
                for i in range(num_c)
            },
        )

        mesh = from_meshio(
            mesh_from_OrderedDict(
                polygons, resolutions, filename="mesh.msh", default_resolution_max=0.4
            )
        )
        # mesh.draw().show()

        basis0 = Basis(mesh, ElementTriP0(), intorder=5)

        epsilon = basis0.zeros() + 3.4777**2
        epsilon[basis0.get_dofs(elements=("box"))] = 1.444**2
        epsilon[basis0.get_dofs(elements=("clad"))] = 1**2
        # basis0.plot(epsilon, colorbar=True).show()
        modes_both = compute_modes(
            basis0, epsilon, wavelength=wavelength, num_modes=2, order=2
        ).sorted(lambda x: -x.n_eff)
        # modes_both[0].show("E", direction="x")
        # modes_both[1].show("E", direction="x")
        coupling_length = wavelength / (2 * np.real(modes_both[0].n_eff - modes_both[1].n_eff))
        # print(f"Maximum power transfer after {coupling_length} um prop length")
        coupling_lengths.append(coupling_length)

    plt.plot(*reference_data, label="Reference")
    plt.plot(w_core_cs * 1e3, coupling_lengths, "ro", label="Calculated")
    plt.legend()
    plt.show()
