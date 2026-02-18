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
# # Modes of waveguides with diagonal anisotropy
# Example based on {cite}`Wang2020`.

# %%
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0, speed_of_light
from shapely import Polygon, box, difference
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

# %% [markdown]
# ## Sweep over propagation angle phi (Fig. 1 reproduction)
# Now that we observe anisotropic behavior, we compare propagation-direction-dependent effective index with the reference.
# Following {cite}`Wang2020`, we sweep the propagation
# angle phi from 0° (Y-propagation) to 90° (Z-propagation) and compute the effective
# indices of the two fundamental modes at 1.31 μm and 1.55 μm wavelengths.
#
# For X-cut LiNbO3, as phi varies:
# - phi=0°: quasi-TE sees n_e (along Z axis), quasi-TM sees n_o
# - phi=90°: both modes see n_o (propagating along Z axis)
#
# The rotation is achieved by transforming the dielectric tensor in the xy-plane.
#
# Note: for propagation angles that are not 0° or 90°, there will also be a non-zero `epsilon_xz` which is not implemented here (it results in a quadratic eigenvalue problem).
# The error in the effective index is not very large for this example, but it could become important with larger birefringence.


# %%
def get_rotated_epsilon_ln(phi_deg, no, ne):
    """
    Get dielectric tensor for X-cut LiNbO3 rotated by angle phi in the xz-plane.

    For X-cut wafer with optical axis (Z crystal) in the chip plane:
    - phi=0°: Y-propagation, extraordinary axis along x (simulation axis)
    - phi=90°: Z-propagation, both in-plane directions see ordinary index

    The dielectric tensor in the waveguide frame (x,y,z) where:
    - x is horizontal (in-plane, perpendicular to propagation)
    - y is vertical (out-of-plane, X crystal axis)
    - z is propagation direction
    """
    phi = np.radians(phi_deg)

    # For X-cut: X crystal axis is vertical (y in simulation)
    # Y and Z crystal axes are in the xz plane
    # At phi=0: z||Y, x||Z -> epsilon_x = n_e^2, epsilon_z = n_o^2
    # At phi=90: z||Z, x||Y -> epsilon_x = n_o^2, epsilon_z = n_e^2

    # The epsilon_y always sees n_o (X crystal direction)
    epsilon_y = no**2

    # Rotation in xz plane: the extraordinary axis rotates
    # epsilon_x and epsilon_z mix according to rotation matrix
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    # At phi=0: epsilon_x = ne^2, epsilon_z = no^2
    # At phi=90: epsilon_x = no^2, epsilon_z = ne^2
    epsilon_x = ne**2 * cos_phi**2 + no**2 * sin_phi**2
    epsilon_z = no**2 * cos_phi**2 + ne**2 * sin_phi**2

    return np.array([epsilon_x, epsilon_y, epsilon_z])


# Sweep parameters
phi_angles = np.linspace(0, 90, 19)  # 0 to 90 degrees in 5-degree steps
wavelengths = [1.31, 1.55]  # micrometers

# Refractive indices at the two wavelengths, using values close to those in the paper
n_values = {
    1.31: {"no": 2.2226, "ne": 2.1458},
    1.55: {"no": 2.2111, "ne": 2.1376},
}

# Store results
results = {wl: {"phi": [], "neff_1": [], "neff_2": []} for wl in wavelengths}

# Run the sweep
print("\n" + "=" * 60)
print("Phi sweep for X-cut LiNbO3 ridge waveguide")
print("=" * 60)

for wavelength in wavelengths:
    no = n_values[wavelength]["no"]
    ne = n_values[wavelength]["ne"]
    print(f"\nWavelength: {wavelength} μm (no={no:.4f}, ne={ne:.4f})")
    print("-" * 50)

    for phi_deg in phi_angles:
        # Get rotated epsilon for LN
        eps_ln = get_rotated_epsilon_ln(phi_deg, no, ne)

        # Build epsilon field
        epsilons_sweep = {
            "core": np.sqrt(eps_ln),  # store as n, will be squared below
            "slab": np.sqrt(eps_ln),
            "box": 1.44,
            "clad": 1.44,
        }

        epsilon_sweep = basis0.zeros()
        for subdomain, e in epsilons_sweep.items():
            epsilon_sweep[basis0.get_dofs(elements=subdomain).all().reshape(-1, 3)] = (
                e**2 if np.isscalar(e) else e**2
            )

        # Compute modes
        modes_sweep = compute_modes(
            basis0, epsilon_sweep, wavelength=wavelength, num_modes=2, order=2
        )

        # Extract effective indices (sorted by real part, descending)
        neffs = sorted([mode.n_eff.real for mode in modes_sweep], reverse=True)

        results[wavelength]["phi"].append(phi_deg)
        results[wavelength]["neff_1"].append(neffs[0])
        results[wavelength]["neff_2"].append(neffs[1])

        print(f"  phi={phi_deg:5.1f}°: n_eff,1 = {neffs[0]:.5f}, n_eff,2 = {neffs[1]:.5f}")

# %% [markdown]
# ### Plot results (reproducing Fig. 1d from the paper)

# %%
# Load reference data from CSV files (extracted from IEEE paper Fig. 1d)
import os

ref_data_dir = os.path.join(os.path.dirname(__file__), "anisotropy")

# Reference data mapping:
# one.csv: mode 1 at 1.31 μm (higher n_eff, solid line in paper)
# two.csv: mode 2 at 1.31 μm (lower n_eff, solid line in paper)
# three.csv: mode 1 at 1.55 μm (higher n_eff, dashed line in paper)
# four.csv: mode 2 at 1.55 μm (lower n_eff, dashed line in paper)

ref_files = {
    "1.31_1": os.path.join(ref_data_dir, "one.csv"),
    "1.31_2": os.path.join(ref_data_dir, "two.csv"),
    "1.55_1": os.path.join(ref_data_dir, "three.csv"),
    "1.55_2": os.path.join(ref_data_dir, "four.csv"),
}

ref_data = {}
for key, filepath in ref_files.items():
    if os.path.exists(filepath):
        data = np.loadtxt(filepath, delimiter=",")
        # Sort by phi angle for proper plotting
        sorted_indices = np.argsort(data[:, 0])
        ref_data[key] = {
            "phi": data[sorted_indices, 0],
            "neff": data[sorted_indices, 1],
        }
        print(f"Loaded reference data: {key} ({len(data)} points)")
    else:
        print(f"Warning: Reference file not found: {filepath}")

# %%
fig, ax = plt.subplots(figsize=(10, 7))

# Plot simulation results
colors = {1.31: ("tab:orange", "tab:blue"), 1.55: ("tab:green", "tab:red")}
linestyles = {1.31: "-", 1.55: "--"}
labels = {1.31: "1.31 μm", 1.55: "1.55 μm"}

for wavelength in wavelengths:
    phi = results[wavelength]["phi"]
    neff_1 = results[wavelength]["neff_1"]
    neff_2 = results[wavelength]["neff_2"]

    ax.plot(
        phi,
        neff_1,
        linestyle=linestyles[wavelength],
        color=colors[wavelength][0],
        label=f"Sim: Mode 1 ({labels[wavelength]})",
        linewidth=2,
    )
    ax.plot(
        phi,
        neff_2,
        linestyle=linestyles[wavelength],
        color=colors[wavelength][1],
        label=f"Sim: Mode 2 ({labels[wavelength]})",
        linewidth=2,
    )

# Overlay reference data with markers
ref_colors = {
    "1.31_1": "tab:orange",
    "1.31_2": "tab:blue",
    "1.55_1": "tab:green",
    "1.55_2": "tab:red",
}
ref_markers = {
    "1.31_1": "o",
    "1.31_2": "s",
    "1.55_1": "^",
    "1.55_2": "v",
}
ref_labels = {
    "1.31_1": "Ref: Mode 1 (1.31 μm)",
    "1.31_2": "Ref: Mode 2 (1.31 μm)",
    "1.55_1": "Ref: Mode 1 (1.55 μm)",
    "1.55_2": "Ref: Mode 2 (1.55 μm)",
}

for key, data in ref_data.items():
    ax.scatter(
        data["phi"],
        data["neff"],
        color=ref_colors[key],
        marker=ref_markers[key],
        s=50,
        label=ref_labels[key],
        edgecolors="none",
    )

ax.set_xlabel("Propagation angle φ (degrees)", fontsize=12)
ax.set_ylabel("Effective refractive index", fontsize=12)
ax.set_title(
    "Effective indices vs propagation angle for X-cut LiNbO3 ridge waveguide\n"
    f"(w={wg_width} μm, t_total={wg_thickness} μm, r_slab={slab_thickness/wg_thickness:.2f})"
)
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 90)

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
