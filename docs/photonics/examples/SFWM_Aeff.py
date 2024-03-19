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
# # Optimisation of spontaneous four-wave mixing in a ring microcavity

# Here we reproduce {cite}`Chuprina2017`

# %% tags=["remove-stderr", "hide-input", "thebe-init"]
from collections import OrderedDict

import numpy as np
from scipy.constants import c, epsilon_0
from shapely.geometry import box
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0, Functional
from skfem.helpers import dot
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

# %% [markdown]
#
# $$
# \int\int\left|u(x,y)\right|^{2}\mathrm{d}x\mathrm{d}y=1
# $$
#
# $$
# A_{\mathrm{eff}}
# =
# \frac{1}{
#   \displaystyle
#   \int\int
#   \mathrm{d}x \mathrm{d}y
#   u_{\mathrm{p}}(x,y)
#   u_{\mathrm{p}}(x,y)
#   u_{\mathrm{s}}^{*}(x,y)
#   u_{\mathrm{i}}^{*}(x,y)
# }
# $$


# %%
def calculate_sfwm_Aeff(basis: Basis, mode_p, mode_s, mode_i) -> np.complex64:
    """
    Calculates the overlap integral for SFWM process by considering the interactions
    between pump, signal, and idler modes in the xy plane.

    Args:
        basis (Basis): Common basis used for all modes in the same geometric structure.
        mode_p, mode_s, mode_i: Mode instances for pump, signal, and idler, respectively.

    Returns:
        np.complex64: The Aeff result for the SFWM process(1/overlap integral).
    """

    def normalization_factor_mode(mode):
        @Functional
        def E2(w):
            return dot(w["E"][0], np.conj(w["E"][0]))

        E = mode.basis.interpolate(mode.E)  # [0]=xy [1]=z
        E_squared_integral = E2.assemble(mode.basis, E=E)
        normalization_factor = 1 / np.sqrt(E_squared_integral)
        return normalization_factor  # Return the normalization factor instead of modifying the mode

    # Apply normalization factors to the electric fields for the overlap calculation
    @Functional(dtype=np.complex64)
    def sfwm_overlap(w):

        E_p_xy = w["E_p"][0]  # shape: (2, x, 3)
        E_s_xy = w["E_s"][0]  # shape: (2, x, 3)
        E_i_xy = w["E_i"][0]  # shape: (2, x, 3)

        overlap_Ex = E_p_xy[0, :, 0] * E_p_xy[0, :, 0] * np.conj(E_s_xy[0, :, 0]) * np.conj(
            E_i_xy[0, :, 0]
        ) + E_p_xy[1, :, 0] * E_p_xy[1, :, 0] * np.conj(E_s_xy[1, :, 0]) * np.conj(E_i_xy[1, :, 0])
        overlap_Ey = E_p_xy[0, :, 1] * E_p_xy[0, :, 1] * np.conj(E_s_xy[0, :, 1]) * np.conj(
            E_i_xy[0, :, 1]
        ) + E_p_xy[1, :, 1] * E_p_xy[1, :, 1] * np.conj(E_s_xy[1, :, 1]) * np.conj(E_i_xy[1, :, 1])
        overlap_Ez = E_p_xy[0, :, 2] * E_p_xy[0, :, 2] * np.conj(E_s_xy[0, :, 2]) * np.conj(
            E_i_xy[0, :, 2]
        ) + E_p_xy[1, :, 2] * E_p_xy[1, :, 2] * np.conj(E_s_xy[1, :, 2]) * np.conj(E_i_xy[1, :, 2])

        return np.array([overlap_Ex, overlap_Ey, overlap_Ez]).T

        # return dot(w["E_p"][0], w["E_p"][0]) * dot(np.conj(w["E_s"][0]), np.conj(w["E_i"][0]))#?
        # return dot(w["E_p"][0], np.conj(w["E_s"][0])) * dot(w["E_p"][0], np.conj(w["E_i"][0]))#??

    overlap_result = sfwm_overlap.assemble(
        basis,
        E_p=mode_p.basis.interpolate(mode_p.E * normalization_factor_mode(mode_p)),
        E_s=mode_s.basis.interpolate(mode_s.E * normalization_factor_mode(mode_s)),
        E_i=mode_i.basis.interpolate(mode_i.E * normalization_factor_mode(mode_i)),
    )

    return 1 / overlap_result


# %% [markdown]
#
# Dispersion relations of materials
# Silicon nitride {cite}`Luke2015`


# %%
def n_X(wavelength):
    x = wavelength
    return (
        1
        + 2.19244563 / (1 - (0.20746607 / x) ** 2)
        + 0.13435116 / (1 - (0.3985835 / x) ** 2)
        + 2.20997784 / (1 - (0.20747044 / x) ** 2)
    ) ** 0.5


# %% [markdown]
#
# Box {cite}`Malitson1965`


# %%
def n_silicon_dioxide(wavelength):
    x = wavelength
    return (
        1
        + 0.6961663 / (1 - (0.0684043 / x) ** 2)
        + 0.4079426 / (1 - (0.1162414 / x) ** 2)
        + 0.8974794 / (1 - (9.896161 / x) ** 2)
    ) ** 0.5


Clad = 1

# %% [markdown]
#
# Waveguide dimensions 1050x500nm

# %%

core = box(0, 0, 1.05, 0.5)
# core = box(0, 0, .5, 0.39)  # 500x390nm
polygons = OrderedDict(
    core=core,
    box=clip_by_rect(core.buffer(1.5, resolution=4), -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(core.buffer(1.5, resolution=4), -np.inf, 0, np.inf, np.inf),
)

resolutions = {"core": {"resolution": 0.025, "distance": 2.0}}

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.2))

num_modes = 2

# For SFWM we have energy conservation and momemtum(k) conservation for 2pumps and signal+idler
# lambda_p0 = 1.4
# lambda_s0 = 1.097
# lambda_i0 = 1.686

lambda_p0 = 1.55
lambda_s0 = 1.55
lambda_i0 = 1.55

basis0 = Basis(mesh, ElementTriP0())

epsilon_p = basis0.zeros()
epsilon_s = basis0.zeros()
epsilon_i = basis0.zeros()


for wavelength, epsilon in zip(
    [lambda_p0, lambda_s0, lambda_i0], [epsilon_p, epsilon_s, epsilon_i]
):
    for subdomain, n_func in {
        "core": n_X,
        "box": n_silicon_dioxide,
        "clad": lambda _: Clad,
    }.items():
        n = n_func(wavelength)
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2


modes_p = compute_modes(basis0, epsilon_p, wavelength=lambda_p0, num_modes=num_modes, order=1)
modes_s = compute_modes(basis0, epsilon_s, wavelength=lambda_s0, num_modes=num_modes, order=1)
modes_i = compute_modes(basis0, epsilon_i, wavelength=lambda_i0, num_modes=num_modes, order=1)

# modes_p[0].show(modes_p[0].E.real)

mode_p = max(modes_p, key=lambda mode: mode.te_fraction)
mode_s = max(modes_s, key=lambda mode: mode.te_fraction)
mode_i = max(modes_i, key=lambda mode: mode.te_fraction)

A_eff = np.real(calculate_sfwm_Aeff(basis0, mode_p, mode_s, mode_i))
print("Aeff in um2:", A_eff)

# Calculation for non-linear coef
chi_3 = 5e-21  # m^2/V^2  #7e-20?
lambda_p0_m = lambda_p0 * 1e-6  # to m
n_p0 = np.real(mode_p.n_eff)
A_eff_m2 = A_eff * 1e-12  # to m^2

omega_p0 = 2 * np.pi * c / lambda_p0_m

# %% [markdown]
# $$
# \gamma=\frac{3\chi^{(3)}\omega_{\mathrm{p}}}{4\varepsilon_{0}c^{2}n_{\mathrm{p0}}^{2}A_{\mathrm{eff}}}
# $$

# %%

gamma = (3 * chi_3 * omega_p0) / (4 * epsilon_0 * c**2 * n_p0**2 * A_eff_m2)

print("gamma:", gamma)

# %% [markdown]
# ## Bibliography
#
# ```{bibliography}
# :style: unsrt
# :filter: docname in docnames
# ```
