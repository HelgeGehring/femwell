# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: femwell
#     language: python
#     name: python3
# ---

# # Photodetector field profiles
#
# ```{caution}
# **It seems that in its current formulation, the detector basis can only represent up to 95% of the input mode (see plot below), indicating it is not complete. This does not occur if the germanium index is set purely real.**
# ```
#
# ```{caution}
# **This example ignores the reflection to cladding at the end of the detector. For long enough detectors where most of the light is absorbed, this should be a small effect.**
# ```
#
# ```{caution}
# **In this example, the scipy solver returns diverging high-order modes. The slepc solver does not seem to have this issue.**
# ```
#
# Mode solving can be used to calculate the optical intensity profile along the length of an absorber. This is useful to calculate optical generation terms for semiconductor simulations.

# +
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from tqdm import tqdm

from femwell.maxwell.waveguide import Mode, calculate_hfield, compute_modes
from femwell.mesh import mesh_from_OrderedDict
from femwell.visualization import plot_domains

# -

# Consider a typical germanium-on-silicon vertical photodetector profile, potentially with some sidewalls. We create a single mesh that has both the silicon input width as well as the larger silicon width under the germanium to avoid mesh interpolation errors:

# +
silicon_core_thickness = 0.22
germanium_thickness = 0.5
germanium_sidewall_angle = 10
clad_vertical_offset = 3

silicon_input_width = 1.5
silicon_detector_width = 3
germanium_width = 1
clad_horizontal_offset = 3

silicon_input = shapely.geometry.box(
    -silicon_input_width / 2, -silicon_core_thickness, silicon_input_width / 2, 0
)
silicon_detector = shapely.geometry.box(
    -silicon_detector_width / 2, -silicon_core_thickness, silicon_detector_width / 2, 0
)
germanium = shapely.Polygon(
    (
        (-germanium_width / 2, 0),
        (
            -germanium_width / 2
            + germanium_thickness / np.tan(np.rad2deg(germanium_sidewall_angle)),
            germanium_thickness,
        ),
        (
            germanium_width / 2
            - germanium_thickness / np.tan(np.rad2deg(germanium_sidewall_angle)),
            germanium_thickness,
        ),
        (germanium_width / 2, 0),
    )
)
clad = shapely.geometry.box(
    -silicon_detector_width / 2 - clad_horizontal_offset,
    -silicon_core_thickness - clad_vertical_offset,
    silicon_detector_width / 2 + clad_horizontal_offset,
    germanium_thickness + clad_vertical_offset,
)

polygons = OrderedDict(
    silicon_input=silicon_input, silicon_detector=silicon_detector, germanium=germanium, clad=clad
)

resolutions = dict(
    silicon_input={"resolution": 0.025, "distance": 5},
    silicon_detector={"resolution": 0.025, "distance": 5},
    germanium={"resolution": 0.025, "distance": 5},
)

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.5))
mesh.draw().show()

plot_domains(mesh)
plt.show()
# -

# By tagging the materials appropriately, we can generate the cross-section of the input silicon:

# +
wavelength = 1.55

# Using refractive indices at 1.55 um
si_index = 3.45
ge_index = 4.6530  # + 1j * 0.29800
clad_index = 1.44

# +
basis0 = Basis(mesh, ElementTriP0(), intorder=4)

epsilon_input = basis0.zeros(dtype=complex)

for subdomain, n in {
    "silicon_input": si_index,
    "silicon_detector": clad_index,
    "germanium": clad_index,
    "clad": clad_index,
}.items():
    epsilon_input[basis0.get_dofs(elements=subdomain)] = n**2

fig, axs = plt.subplots(1, 2)
for ax in axs:
    ax.set_aspect(1)
axs[0].set_title(r"$\Re\epsilon$, input")
basis0.plot(epsilon_input.real, colorbar=True, ax=axs[0])
axs[1].set_title(r"$\Im\epsilon$, input")
basis0.plot(epsilon_input.imag, shading="gouraud", colorbar=True, ax=axs[1])
plt.show()
# -

# The fundamental mode of this geometry will be our input mode:

input_modes = compute_modes(
    basis0,
    epsilon_input,
    wavelength=wavelength,
    num_modes=1,
    order=1,
    radius=np.inf,
    solver="slepc",
)

# +
input_mode = input_modes[0]

input_mode.plot(input_mode.E.real, colorbar=True, direction="x")
plt.show()

# +
basis0 = Basis(mesh, ElementTriP0(), intorder=4)

epsilon_detector = basis0.zeros(dtype=complex)

for subdomain, n in {
    "silicon_input": si_index,
    "silicon_detector": si_index,
    "germanium": ge_index,
    "clad": clad_index,
}.items():
    epsilon_detector[basis0.get_dofs(elements=subdomain)] = n**2

fig, axs = plt.subplots(1, 2)
for ax in axs:
    ax.set_aspect(1)
axs[0].set_title(r"$\Re\epsilon$, detector")
basis0.plot(epsilon_detector.real, colorbar=True, ax=axs[0])
axs[1].set_title(r"$\Im\epsilon$, detector")
basis0.plot(epsilon_detector.imag, shading="gouraud", colorbar=True, ax=axs[1])
plt.show()
# -

detector_modes = compute_modes(
    basis0,
    epsilon_detector,
    wavelength=wavelength,
    num_modes=40,
    order=1,
    radius=np.inf,
    solver="slepc",
)

for i, mode in enumerate(detector_modes):
    if not i % 5:
        print(f"Mode index: {i}")
        mode.plot(mode.E.real, colorbar=True, direction="x")
        plt.show()

#

# +
overlaps = []
ks = []

sum_modes = detector_modes[0].basis.zeros(dtype=complex)

for mode in tqdm(detector_modes):
    overlaps.append(mode.calculate_overlap(input_mode))
    ks.append(mode.k)

    sum_modes += mode.E * mode.calculate_overlap(input_mode)
# -

detector_modes[0].show(sum_modes.real, colorbar=True)

# Here we evaluate the completeness of the detector basis, and see that it fails to capture about 5% of the incoming mode. We can also see the effect of modes with different symmetry and polarization to the input mode (plateaus in cumulative overlap).

plt.plot(np.cumsum(np.abs(overlaps) ** 2))
plt.ylim([0, 1.1])
plt.axhline(y=1, color="k", linestyle="--")
plt.xlabel("Mode index")
plt.ylabel("Cumulative power overlap")

# Following eigenmode expansion, the propagating field inside the detector cross-section given a pure input mode can be expressed as (ignoring reflections):
#
# $$ E(x,y,z) = \sum_k \braket{E_k|E_0} e^{-i \beta_k z} \ket{E_k} $$
#
# where $\braket{\bm{x}|E_k} = E_k(x,y)$ are the cross-sectional mode profiles.

# +
z = 0


def field_profile_at_z(z):
    field = detector_modes[0].basis.zeros(dtype=complex)
    for i, (k, overlap, mode) in enumerate(zip(ks, overlaps, detector_modes)):
        field += overlap * np.exp(-1j * (np.real(k) - 1j * np.abs(np.imag(k))) * z) * mode.E
    return field


# -

for z in np.linspace(0, 2, 11):
    field = field_profile_at_z(z)
    detector_modes[0].show(field.real, colorbar=True)

#
