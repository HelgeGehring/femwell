# -*- coding: utf-8 -*-
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
# # Perturbation theory vs. finite differences
#
# `Mode.calculate_pertubated_neff` implements first-order perturbation theory
# (a Hellmann-Feynman-style formula) to estimate how a mode's effective index
# changes for a given permittivity change $\Delta\varepsilon$, from a single
# mode solve — no re-meshing or re-solving required. This is useful whenever
# you need the sensitivity of $n_{eff}$ to a design parameter (e.g. for
# gradient-based optimization of a waveguide geometry).
#
# This example validates that prediction against literal finite-difference
# re-solves for a width perturbation of a strip waveguide, and shows the
# corresponding speed/accuracy trade-off: perturbation theory reuses a single
# baseline mode solve and adds one cheap surface integral, while finite
# differences need a full extra mesh + solve for every perturbation size.
#
# The key subtlety demonstrated here is how to build $\Delta\varepsilon$ for
# a *shape* (boundary-position) perturbation, as opposed to a perturbation
# where the material itself changes everywhere at once: only the thin sliver
# of material that the moving boundary sweeps over should contribute, and
# nowhere else. We mesh that sliver as its own named subdomain right at the
# sidewall being perturbed, sitting outside the original core.

# %% tags=["hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

# %% [markdown]
# ## Baseline mode
#
# Same strip waveguide as in {doc}`vary_width`: a 500 nm x 220 nm silicon
# core on a uniform SiO2 clad. We additionally mesh a thin "sliver"
# subdomain right outside the core's right sidewall -- at the baseline
# width this sliver is physically just more cladding (same index as the
# surrounding clad), but it marks exactly the region that would turn into
# core material if the waveguide were widened by `dw`.

# %% tags=["remove-stderr"]
wavelength = 1.55
width = 0.5
thickness = 0.22
n_core, n_clad = 3.48, 1.444
dw = 0.02

core = box(0, 0, width, thickness)
sliver = box(width, 0, width + dw, thickness)
clad = core.buffer(1.0, resolution=6)

polygons = OrderedDict(sliver=sliver, core=core, clad=clad)
resolutions = {
    "core": {"resolution": 0.04, "distance": 0.5},
    "sliver": {"resolution": dw / 3, "distance": 0.2},
}

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.3))
basis0 = Basis(mesh, ElementTriP0())

epsilon_baseline = basis0.zeros(dtype=complex)
for subdomain, n in {"core": n_core, "sliver": n_clad, "clad": n_clad}.items():
    epsilon_baseline[basis0.get_dofs(elements=subdomain)] = n**2

mode_baseline = compute_modes(basis0, epsilon_baseline, wavelength=wavelength, num_modes=1, order=1)[0]
print(f"baseline n_eff = {mode_baseline.n_eff.real:.6f}")

# %% [markdown]
# ## Perturbation-theory prediction
#
# $\Delta\varepsilon$ is zero everywhere except inside the sliver, where it
# equals $n_{core}^2 - n_{clad}^2$ -- exactly the material change that
# widening the waveguide by `dw` would cause there.

# %% tags=["remove-stderr"]
delta_epsilon = basis0.zeros(dtype=complex)
delta_epsilon[basis0.get_dofs(elements="sliver")] = n_core**2 - n_clad**2

n_eff_theory = mode_baseline.calculate_pertubated_neff(delta_epsilon).real
print(f"perturbation-theory predicted n_eff = {n_eff_theory:.6f}")

# %% [markdown]
# ## Finite-difference ground truth
#
# Mesh and solve the actual, independently-built `width + dw` geometry from
# scratch -- the expensive way, exactly as {doc}`vary_width` does for every
# point of its sweep.

# %% tags=["remove-stderr"]
core_new = box(0, 0, width + dw, thickness)
clad_new = core_new.buffer(1.0, resolution=6)
polygons_new = OrderedDict(core=core_new, clad=clad_new)
mesh_new = from_meshio(
    mesh_from_OrderedDict(polygons_new, {"core": {"resolution": 0.04, "distance": 0.5}}, default_resolution_max=0.3)
)
basis0_new = Basis(mesh_new, ElementTriP0())
epsilon_new = basis0_new.zeros(dtype=complex)
for subdomain, n in {"core": n_core, "clad": n_clad}.items():
    epsilon_new[basis0_new.get_dofs(elements=subdomain)] = n**2

mode_new = compute_modes(basis0_new, epsilon_new, wavelength=wavelength, num_modes=1, order=1)[0]
n_eff_fd = mode_new.n_eff.real
print(f"finite-difference n_eff      = {n_eff_fd:.6f}")

# %% [markdown]
# ## Agreement

# %% tags=["remove-stderr"]
dneff_theory = n_eff_theory - mode_baseline.n_eff.real
dneff_fd = n_eff_fd - mode_baseline.n_eff.real
rel_diff = abs(dneff_theory - dneff_fd) / abs(dneff_fd) * 100

print(f"dn_eff (perturbation theory) = {dneff_theory:.6e}")
print(f"dn_eff (finite difference)   = {dneff_fd:.6e}")
print(f"relative difference          = {rel_diff:.2f}%")

# %% [markdown]
# ## Convergence with perturbation size
#
# Perturbation theory is a first-order (linear) approximation, so it should
# get *more* accurate as the perturbation shrinks, while a plain finite
# difference eventually loses precision to mesh/discretization noise if
# pushed too small. We sweep `dw` to show both trends on the same axes.

# %% tags=["hide-input", "remove-stderr"]
dws = np.array([0.02, 0.01, 0.005])
rel_diffs = np.zeros_like(dws)

for i, dw_i in enumerate(dws):
    sliver_i = box(width, 0, width + dw_i, thickness)
    polygons_i = OrderedDict(sliver=sliver_i, core=core, clad=clad)
    resolutions_i = {
        "core": {"resolution": 0.04, "distance": 0.5},
        "sliver": {"resolution": dw_i / 3, "distance": 0.2},
    }
    mesh_i = from_meshio(mesh_from_OrderedDict(polygons_i, resolutions_i, default_resolution_max=0.3))
    basis0_i = Basis(mesh_i, ElementTriP0())
    epsilon_i = basis0_i.zeros(dtype=complex)
    for subdomain, n in {"core": n_core, "sliver": n_clad, "clad": n_clad}.items():
        epsilon_i[basis0_i.get_dofs(elements=subdomain)] = n**2
    mode_i = compute_modes(basis0_i, epsilon_i, wavelength=wavelength, num_modes=1, order=1)[0]

    delta_epsilon_i = basis0_i.zeros(dtype=complex)
    delta_epsilon_i[basis0_i.get_dofs(elements="sliver")] = n_core**2 - n_clad**2
    n_eff_theory_i = mode_i.calculate_pertubated_neff(delta_epsilon_i).real

    core_new_i = box(0, 0, width + dw_i, thickness)
    clad_new_i = core_new_i.buffer(1.0, resolution=6)
    polygons_new_i = OrderedDict(core=core_new_i, clad=clad_new_i)
    mesh_new_i = from_meshio(
        mesh_from_OrderedDict(
            polygons_new_i, {"core": {"resolution": 0.04, "distance": 0.5}}, default_resolution_max=0.3
        )
    )
    basis0_new_i = Basis(mesh_new_i, ElementTriP0())
    epsilon_new_i = basis0_new_i.zeros(dtype=complex)
    for subdomain, n in {"core": n_core, "clad": n_clad}.items():
        epsilon_new_i[basis0_new_i.get_dofs(elements=subdomain)] = n**2
    mode_new_i = compute_modes(basis0_new_i, epsilon_new_i, wavelength=wavelength, num_modes=1, order=1)[0]

    dneff_theory_i = n_eff_theory_i - mode_i.n_eff.real
    dneff_fd_i = mode_new_i.n_eff.real - mode_i.n_eff.real
    rel_diffs[i] = abs(dneff_theory_i - dneff_fd_i) / abs(dneff_fd_i) * 100
    print(f"dw={dw_i * 1000:6.1f}nm: dn_eff theory={dneff_theory_i:.4e}  fd={dneff_fd_i:.4e}  rel.diff={rel_diffs[i]:.2f}%")

# %% tags=["hide-input"]
plt.figure()
plt.loglog(dws * 1000, rel_diffs, "-o")
plt.xlabel("perturbation size dw / nm")
plt.ylabel("relative difference vs. finite difference / %")
plt.title("Perturbation theory converges to finite differences as dw -> 0")
plt.grid(True, which="both")
plt.show()

# %% [markdown]
# The relative disagreement shrinks as `dw` shrinks, confirming that
# `calculate_pertubated_neff` is genuinely computing the correct first-order
# derivative, not just landing close by coincidence for one arbitrary `dw`.
# In particular, for a gradient-based optimization loop that repeatedly needs
# $dn_{eff}/dw$ at many candidate widths, using perturbation theory to
# extract the derivative directly (one baseline solve + one cheap surface
# integral) is dramatically cheaper than the finite-difference approach used
# by default in {doc}`vary_width` -- no second re-mesh/re-solve is needed at
# all.

# %%
