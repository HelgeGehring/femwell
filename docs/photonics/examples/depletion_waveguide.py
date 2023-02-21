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

# # PN junction modulator

# + tags=["remove-stderr", "hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.affinity
from scipy.constants import epsilon_0, speed_of_light
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver import (
    calculate_hfield,
    calculate_overlap,
    compute_modes,
    plot_mode,
)
from femwell.pn_analytical import *

# -

# We can study the propagation constant in waveguides as a function of arbitrary physics.
# Here, we consider the depletion approximation to pn junctions to study how doping level and junction placement affect modulation in a doped silicon waveguide.

# +
clad_thickness = 2
slab_width = 3
wg_width = 0.5
wg_thickness = 0.22
slab_thickness = 0.09
core = shapely.geometry.box(-wg_width / 2, -wg_thickness / 2, wg_width / 2, wg_thickness / 2)
slab = shapely.geometry.box(
    -slab_width / 2, -wg_thickness / 2, slab_width / 2, -wg_thickness / 2 + slab_thickness
)
clad = shapely.geometry.box(
    -slab_width / 2, -clad_thickness / 2, slab_width / 2, clad_thickness / 2
)

polygons = OrderedDict(
    core=core,
    slab=slab,
    clad=clad,
)

resolutions = dict(
    core={"resolution": 0.02, "distance": 0.5}, slab={"resolution": 0.04, "distance": 0.5}
)

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=10))
mesh.draw().show()
# -

# To define the epsilon, we proceed as for a regular waveguide, but we superimpose a voltage-dependent index of refraction based on the Soref Equations. These phenomenologically relate the change in complex index of refraction of silicon as a function of the concentration of free carriers. We model the spatial distribution of carriers according to the physics of a 1D PN junction in the depletion approximation. For more accurate results, full modeling of the silicon processing and physics through TCAD must be performed.

# # +
xpn = 0
NA = 1e17
ND = 1e17
V = 0
wavelength = 1.55


def define_index(mesh, V, xpn, NA, ND, wavelength):
    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    for subdomain, n in {"core": 3.45, "slab": 3.45}.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2
    epsilon += basis0.project(
        lambda x: index_pn_junction(x[0], xpn, NA, ND, V, wavelength) ** 2,
        dtype=complex,
    )
    for subdomain, n in {"clad": 1.444}.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2
    return basis0, epsilon


basis0, epsilon = define_index(mesh, V, xpn, NA, ND, wavelength)
basis0.plot(epsilon.real, colorbar=True).show()
basis0.plot(epsilon.imag, colorbar=True).show()

# -
# The index change is weak compared to the contrast between silicon and silicon dioxide, but it is accompanied by a change in absorption which is easier to observe. As voltage is increased, the region without charge widens, which is the mechanism behind depletion modulation:

V = -4
basis0, epsilon = define_index(mesh, V, xpn, NA, ND, wavelength)
basis0.plot(epsilon.imag, colorbar=True).show()
# -

# And now we can mode solve as before, and observe the change in effective index and absorption of the fundamental mode as a function of applied voltage for given junction position and doping levels:

# +
voltages = [0, -1, -2, -3, -4]
neff_vs_V = []
for V in voltages:
    basis0, epsilon = define_index(mesh, V, xpn, NA, ND, wavelength)
    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=1, order=2)
    neff_vs_V.append(lams)

plt.plot(voltages, neff_vs_V.real)
plt.title(f"NA = {NA}, ND = {ND}, xpn = {xpn}, wavelength = {wavelength}")
plt.xlabel("Voltage (V)")
plt.ylabel("neff0")

plt.plot(voltages, k_to_alpha_dB(neff_vs_V.imag, wavelength))
plt.title(f"NA = {NA}, ND = {ND}, xpn = {xpn}, wavelength = {wavelength}")
plt.xlabel("Voltage (V)")
plt.ylabel("absorption (dB/cm)")

# -

# References:

# From Chrostowski, L., & Hochberg, M. (2015). Silicon Photonics Design: From Devices to Systems. Cambridge University Press. doi: 10.1017/CBO9781316084168
#     Citing:
#     (1) R. Soref and B. Bennett, "Electrooptical effects in silicon," in IEEE Journal of Quantum Electronics, vol. 23, no. 1, pp. 123-129, January 1987, doi: 10.1109/JQE.1987.1073206.
#     (2) Reed, G. T., Mashanovich, G., Gardes, F. Y., & Thomson, D. J. (2010). Silicon optical modulators. Nature Photonics, 4(8), 518–526. doi: 10.1038/nphoton.2010.179
#     (3) M. Nedeljkovic, R. Soref and G. Z. Mashanovich, "Free-Carrier Electrorefraction and Electroabsorption Modulation Predictions for Silicon Over the 1–14- $\mu\hbox{m}$ Infrared Wavelength Range," in IEEE Photonics Journal, vol. 3, no. 6, pp. 1171-1180, Dec. 2011, doi: 10.1109/JPHOT.2011.2171930.
