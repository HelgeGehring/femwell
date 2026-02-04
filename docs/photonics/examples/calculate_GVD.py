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
# # Calculate GVD of waveguide

# Reproduce figure 2e from {cite}`Klenner2016`


# %% tags=["remove-stderr", "hide-input", "thebe-init"]
import math
import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from meshwell.cad import cad
from meshwell.mesh import mesh
from meshwell.polysurface import PolySurface
from meshwell.resolution import ThresholdField
from scipy.interpolate import UnivariateSpline
from skfem import Basis, ElementTriP0
from skfem.io import from_meshio
from tqdm import tqdm

from femwell.maxwell.waveguide import compute_modes
from femwell.visualization import plot_domains

# %% [markdown]
# First, construct the geometry of waveguide according to figure 2d of the paper.
# %%
width = 0.88  # um
height = 0.69  # um

wavelength_range = [500, 2500]
wavelegnth_step = 50

# Construct waveguide geometry
core = shapely.geometry.box(-width / 2, 0, +width / 2, height)
cladding = shapely.geometry.box(-width * 2, 0, width * 2, height * 3)
buried_oxide = shapely.geometry.box(-width * 2, -height * 2, width * 2, 0)

# Create meshwell PolySurface entities with correct mesh_order (core=0, cladding=1, buried_oxide=2)
core_surface = PolySurface(polygons=core, physical_name="core", mesh_order=0)
cladding_surface = PolySurface(polygons=cladding, physical_name="cladding", mesh_order=1)
buried_oxide_surface = PolySurface(
    polygons=buried_oxide, physical_name="buried_oxide", mesh_order=2
)

# Define resolution specifications using ThresholdField
resolution_specs = {
    "core": [ThresholdField(sizemin=0.02, distmax=20, sizemax=0.5, apply_to="surfaces")],
}

# Create entities list
entities = [core_surface, buried_oxide_surface, cladding_surface]

# Generate CAD and mesh, using meshwell functional interface
with tempfile.TemporaryDirectory() as tmp:
    cad_path = os.path.join(tmp, "waveguide.xao")
    msh_path = os.path.join(tmp, "waveguide.msh")
    cad(entities_list=entities, output_file=cad_path)
    output_mesh = mesh(
        dim=2,
        input_file=cad_path,
        output_file=msh_path,
        default_characteristic_length=0.5,
        resolution_specs=resolution_specs,
        n_threads=1,
    )

# %%
# Load mesh for femwell
mesh = from_meshio(output_mesh)
mesh.draw().show()
plot_domains(mesh)
plt.show()


# %% [markdown]
# Use sellmeier equation to determine the refractive index of material.


# %%
def n_Si3N4(wavelength, fit=False):
    if wavelength >= 0.31 and wavelength <= 5.507:
        return math.sqrt(
            (3.0249 * wavelength**2) / (wavelength**2 - 0.1353406**2)
            + (40314 * wavelength**2) / (wavelength**2 - 1239.842**2)
            + 1
        )
    else:
        raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for Si3N4")


def n_SiO2(wavelength):
    if wavelength < 0.21 or wavelength > 6.7:
        raise ValueError(f"wavelength provided is {wavelength}um, is out of the range for {type}")
    return np.sqrt(
        0.6961663 * wavelength**2 / (wavelength**2 - 0.0684043**2)
        + (0.4079426 * wavelength**2 / (wavelength**2 - 0.1162414**2))
        + (0.8974794 * wavelength**2 / (wavelength**2 - 9.896161**2))
        + 1
    )


n_dict = {"core": n_Si3N4, "cladding": n_SiO2, "buried_oxide": n_SiO2}

# %% [markdown]
# Create the mesh, and sweep wavelength using the same mesh. The target mode is te mode, so the mode is selected by highest te fraction
# %% tags=["remove-stderr"]
neff_list = []
aeff_list = []
basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros()
wavelength_list = np.linspace(wavelength_range[0], wavelength_range[1], wavelegnth_step)

for wavelength in tqdm(wavelength_list):
    wavelength = wavelength * 1e-3
    for subdomain, n in n_dict.items():
        epsilon[basis0.get_dofs(elements=subdomain)] = n(wavelength) ** 2

    modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=3, order=1)
    modes_sorted = modes.sorted(key=lambda mode: -np.real(mode.te_fraction))
    mode = modes_sorted[0]

    neff_list.append(np.real(mode.n_eff))
    aeff_list.append(np.real(mode.calculate_effective_area()))

# %% [markdown]
# Calculate the GVD by fitting a curve for wavelength vs neff. Then take second derivative of the curve
# %%

y_spl = UnivariateSpline(wavelength_list, neff_list, s=0, k=3)
x_range = np.linspace(wavelength_list[0], wavelength_list[-1], 1000)
y_spl_2d = y_spl.derivative(n=2)

# Plot the result
fig, axs = plt.subplots(3, 1, figsize=(9, 20))

axs[0].set_xlabel("Wavelength / nm")
axs[0].set_ylabel("neff")
axs[0].set_title(" neff vs wavelength fit")
axs[0].semilogy(x_range, y_spl(x_range))
axs[0].semilogy(wavelength_list, neff_list, "ro", label="data")
axs[0].legend()
axs[0].set_xlim(500, 2200)

axs[1].set_xlabel("Wavelength / nm")
axs[1].set_ylabel("neff''")
axs[1].set_title("wavelength vs second derivative of neff")
axs[1].plot(x_range, y_spl_2d(x_range))
axs[1].set_xlim(500, 2200)

# ----plot reference data-------
ref_gvd = pd.read_csv("../reference_data/Klenner/GVD.csv", dtype=np.float64)
ref_gvd_x, ref_gvd_y = np.split(ref_gvd.values, 2, axis=1)
axs[2].plot(ref_gvd_x, ref_gvd_y, c="green", label="paper")

# ----Calculate and plot GVD
GVD = -wavelength_list / (2.99792e-7) * y_spl_2d(wavelength_list)
axs[2].scatter(wavelength_list, GVD, label="calculated", c="red")

axs[2].set_ylabel("GVD")
axs[2].set_xlabel("Wavelength / nm")
axs[2].set_ylim(-1000, 200)
axs[2].set_xlim(500, 2200)
axs[2].set_title("GVD paramter")
axs[2].legend()

plt.tight_layout()
plt.show()

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```

# %%
