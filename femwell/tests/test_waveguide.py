from collections import OrderedDict

import numpy as np
import pytest
from shapely.geometry import box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict


def test_calculate_pertubated_neff_matches_finite_difference():
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

    delta_epsilon = basis0.zeros(dtype=complex)
    delta_epsilon[basis0.get_dofs(elements="sliver")] = n_core**2 - n_clad**2
    n_eff_theory = mode_baseline.calculate_pertubated_neff(delta_epsilon).real
    dneff_theory = n_eff_theory - mode_baseline.n_eff.real

    core_new = box(0, 0, width + dw, thickness)
    clad_new = core_new.buffer(1.0, resolution=6)
    polygons_new = OrderedDict(core=core_new, clad=clad_new)
    mesh_new = from_meshio(
        mesh_from_OrderedDict(
            polygons_new, {"core": {"resolution": 0.04, "distance": 0.5}}, default_resolution_max=0.3
        )
    )
    basis0_new = Basis(mesh_new, ElementTriP0())
    epsilon_new = basis0_new.zeros(dtype=complex)
    for subdomain, n in {"core": n_core, "clad": n_clad}.items():
        epsilon_new[basis0_new.get_dofs(elements=subdomain)] = n**2
    mode_new = compute_modes(basis0_new, epsilon_new, wavelength=wavelength, num_modes=1, order=1)[0]
    dneff_fd = mode_new.n_eff.real - mode_baseline.n_eff.real

    assert dneff_theory == pytest.approx(dneff_fd, rel=0.2)
