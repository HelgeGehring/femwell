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
# # Microstrip waveguide - vary width

# In this example we calculate effective epsilon of the microstrip waveguides from {cite}`Jansen1978`

# %% tags=["hide-input"]

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import shapely
import shapely.ops
from shapely.geometry import LineString, box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from tqdm import tqdm

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict


# %%
def mesh_waveguide(filename, wsim, hclad, hsi, wcore, hcore):
    core = box(-wcore / 2, -hcore / 2, wcore / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)

    combined = shapely.ops.unary_union([clad, silicon])

    polygons = OrderedDict(
        surface=LineString(combined.exterior),
        core_interface=core.buffer(hcore / 10).exterior,
        core=core,
        clad=clad,
        silicon=silicon,
    )

    resolutions = dict(
        core={"resolution": min(wcore / 20, hcore / 4), "distance": 1, "SizeMax": 0.5},
        core_interface={
            "resolution": min(wcore / 20, hcore / 4),
            "distance": 1,
            "SizeMax": 0.5,
        },
        # silicon={"resolution": 0.2, "distance": 5},
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=5)


# %% tags=["hide-output"]
frequencies = np.linspace(1e9, 16e9, 16)
# widths = (0.04, 0.6, 0.8, 0.1, 0.14, 0.18, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.6, 3,)
widths = (
    0.04,
    0.4,
    3,
)
epsilon_effs = np.zeros((len(widths), len(frequencies), 1), dtype=complex)
characteristic_impedances = np.zeros((len(widths), len(frequencies), 1), dtype=complex)

for i, width in enumerate(tqdm(widths)):
    mesh = from_meshio(
        mesh_waveguide(
            filename="mesh.msh",
            wsim=10,
            hclad=10,
            hsi=0.64,
            wcore=width,
            hcore=0.005,
        )
    )
    mesh = mesh.scaled((1e-3,) * 2)
    core_interface = box(-width / 2, -0.005 / 2, width / 2, 0.005 / 2).buffer(0.005 / 10)

    for j, frequency in enumerate(tqdm(frequencies, leave=False)):
        basis0 = Basis(mesh, ElementTriP0(), intorder=4)
        epsilon = basis0.zeros().astype(complex)
        epsilon[basis0.get_dofs(elements="silicon")] = 9.9 + 0.0005
        epsilon[basis0.get_dofs(elements="clad")] = 1.0
        epsilon[basis0.get_dofs(elements="core")] = 1 - 1j * 1 / (
            18e-6 * 1e-3
        ) / scipy.constants.epsilon_0 / (2 * np.pi * frequency)
        # basis0.plot(np.real(epsilon), colorbar=True).show()

        modes = compute_modes(
            basis0,
            epsilon,
            wavelength=scipy.constants.speed_of_light / frequency,
            num_modes=1,
            metallic_boundaries=True,
        )
        print(f"Width: {width}, Frequency: {frequency/1e9} GHz")
        print(f"Effective epsilons {modes.n_effs**2}")
        # modes[0].show("E", part="real", plot_vectors=True, colorbar=True)

        epsilon_effs[i, j] = modes.n_effs**2

        # Power-current impedance: the current is the circulation of H on a contour just outside
        # of the strip (the mesh coordinates are in m, the shapely contour is in mm)
        characteristic_impedances[i, j] = modes[0].calculate_impedance(
            "PI",
            elements=lambda x: shapely.contains_xy(core_interface, x[0] * 1e3, x[1] * 1e3),
        )
        print("characteristic impedance", characteristic_impedances[i, j])

# %% [markdown]
# The characteristic impedance is the power-current impedance $Z_{PI} = 2P/|I|^2$.
# The current $I$ follows from Ampère's law on a contour just outside of the strip.
# The other definitions are available from the same mode, e.g.
# `modes[0].calculate_impedance("PV", start=(0, 0), end=(0, -(0.005 / 2 + 0.64) * 1e-3))`
# for the power-voltage impedance, and `"VI"` (with `elements=`) for the voltage-current
# impedance. For quasi-TEM lines like the microstrip they differ from each other; the
# Jansen-Kirschning model of the microstrip characteristic impedance uses the power-current form.

# %% tags=["hide-input"]
plt.figure(figsize=(10, 14))
plt.xlabel("Frequency / GHz")
plt.ylabel("Effective dielectric constant")

for i, width in enumerate(widths):
    plt.plot(frequencies / 1e9, epsilon_effs[i, :, 0].real)
    plt.annotate(
        xy=(frequencies[-1] / 1e9, epsilon_effs[i, :, 0].real[-1]),
        text=str(width),
        va="center",
    )

plt.show()

plt.figure(figsize=(10, 14))
plt.xlabel("Frequency / GHz")
plt.ylabel("Characteristic Impedance / Ohm")

for i, width in enumerate(widths):
    plt.plot(frequencies / 1e9, characteristic_impedances[i, :, 0].real)
    plt.annotate(
        xy=(frequencies[-1] / 1e9, characteristic_impedances[i, :, 0].real[-1]),
        text=str(width),
        va="center",
    )

plt.show()
# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
