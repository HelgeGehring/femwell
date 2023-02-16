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

# # Modes of a bent rectangular waveguide (in work)

# + tags=["hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import shapely
import shapely.affinity
from shapely.ops import clip_by_rect
from skfem import Basis, ElementDG, ElementTriP1
from skfem.io.meshio import from_meshio
from tqdm import tqdm

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver import (
    calculate_hfield,
    calculate_overlap,
    compute_modes,
    plot_mode,
)

# -

# We describe the geometry using shapely.
# In this case it's simple: we use a shapely.box for the waveguide.
# For the surrounding we buffer the core and clip it to the part below the waveguide for the box.
# The remaining buffer is used as the clad.
# For the core we set the resolution to 30nm and let it fall of over 500nm

# +
wavelength = 1.55

wg_width = 0.5
wg_thickness = 0.22
slab_thickness = 0.11
core = shapely.geometry.box(-wg_width / 2, 0, +wg_width / 2, wg_thickness)
slab = shapely.geometry.box(-1, 0, 4, slab_thickness)
# env = shapely.affinity.translate(core.buffer(4.5, resolution=8), 1.5, -1.5)
env = shapely.geometry.box(-1, -4, 4, wg_thickness + 2)

polygons = OrderedDict(
    core=core,
    slab=slab,
    box=clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
)

resolutions = dict(
    core={"resolution": 0.03, "distance": 1}, slab={"resolution": 0.1, "distance": 0.5}
)

mesh = from_meshio(
    mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.2, filename="mesh.msh")
)
mesh.draw().show()
# -

# On this mesh, we define the epsilon. We do this by setting domainwise the epsilon to the squared refractive index.
# We additionally add a PML layer bt adding a imaginary part to the epsilon

# +
basis0 = Basis(mesh, ElementDG(ElementTriP1()))
epsilon = basis0.zeros(dtype=complex)
for subdomain, n in {"core": 3.48, "slab": 3.48, "box": 1.48, "clad": 1.0}.items():
    epsilon[basis0.get_dofs(elements=subdomain)] = n**2
epsilon += basis0.project(
    lambda x: -10j * (np.maximum(0, x[0] - 2.0) ** 2 + np.maximum(0, -x[1] - 2.0) ** 2),
    dtype=complex,
)
basis0.plot(epsilon.real, shading="gouraud", colorbar=True).show()
basis0.plot(epsilon.imag, shading="gouraud", colorbar=True).show()
# -

# We calculate now the modes for the geometry we just set up.
# We do it for the case, where the bend-radius is infinite, i.e. a straight waveguide
# and for waveguides with a smaller radius.
# Subsequently, we calculate the overlap integrals between the modes to determine the coupling efficiency
# And determine from the imaginary part the bend loss

# + tags=["remove-stderr"]
lams_straight, basis_straight, xs_straight = compute_modes(
    basis0, epsilon, wavelength=wavelength, mu_r=1, num_modes=1, order=2, radius=np.inf
)
H_straight = calculate_hfield(
    basis_straight,
    xs_straight[0],
    2 * np.pi / wavelength * lams_straight[0],
    omega=2 * np.pi / wavelength * scipy.constants.speed_of_light,
)

radiuss = np.linspace(25, 5, 21)
radiuss_lams = []
overlaps = []
lam_guess = lams_straight[0]
for radius in tqdm(radiuss):
    lams, basis, xs = compute_modes(
        basis0,
        epsilon,
        wavelength=wavelength,
        mu_r=1,
        num_modes=1,
        order=2,
        radius=radius,
        n_guess=lam_guess,
    )
    lam_guess = lams[0]
    H_bent = calculate_hfield(
        basis_straight,
        xs[0],
        2 * np.pi / wavelength * lams[0],
        omega=2 * np.pi / wavelength * scipy.constants.speed_of_light,
    )
    radiuss_lams.append(lams[0])

    overlaps.append(
        calculate_overlap(basis_straight, xs[0], H_bent, basis_straight, xs_straight[0], H_straight)
    )

plt.xlabel("Radius / μm")
plt.ylabel("Mode overlap with straight waveguide mode / dB")
plt.plot(radiuss, -10 * np.log10(np.abs(overlaps) ** 2))
plt.show()
plt.xlabel("Radius / μm")
plt.ylabel("Logarithm of imaginary part of refractive index")
plt.plot(radiuss, np.log10(np.abs(np.imag(radiuss_lams))))
plt.show()
plt.xlabel("Radius / μm")
plt.ylabel("90-degree bend transmission / dB")
plt.plot(
    radiuss,
    -10
    * np.log10(
        np.exp(-2 * np.pi / wavelength * radius * np.abs(np.imag(radiuss_lams)) * np.pi / 2)
    ),
)
plt.show()
# -

# We now plot the mode calculated for the smallest bend radius to check that it's still within the waveguide.
# As modes can have complex fields as soon as the epsilon gets complex, so we get a complex field for each mode.
# Here we show only the real part of the mode.

# +
for i, lam in enumerate(lams):
    print(f"Effective refractive index: {lam:.14f}")
    plot_mode(basis, xs[i].real, colorbar=True, direction="x")
    plt.show()
