import tempfile
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import shapely
import shapely.ops
from shapely.geometry import LineString, box
from skfem import Basis, ElementTriP0, Mesh

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict


def mesh_waveguide(filename, wsim, hclad, hsi, wcore, hcore):
    core = box(-wcore / 2, -hcore / 2, wcore / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)

    combined = shapely.ops.unary_union([core, clad, silicon])

    polygons = OrderedDict(
        boundary=LineString(combined.exterior),
        interface=LineString(
            [
                (-wsim / 2, -hcore / 2),
                (wsim / 2, -hcore / 2),
            ]
        ),
        core_interface=core.exterior,
        core=core,
        clad=clad,
        silicon=silicon,
    )

    resolutions = dict(
        core_interface={"resolution": 0.003, "distance": 1},
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=1)


def mesh_waveguide_1(filename, wsim, hclad, hsi, wcore, hcore, gap):
    core_l = box(-wcore - gap, -hcore / 2, -gap / 2, hcore / 2)
    core_r = box(wcore + gap, -hcore / 2, gap / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)

    combined = shapely.ops.unary_union([core_l, core_r, clad, silicon])

    polygons = OrderedDict(
        surface=LineString(combined.exterior),
        interface=LineString(
            [
                (-wsim / 2, -hcore / 2),
                (wsim / 2, -hcore / 2),
            ]
        ),
        core_l_interface=core_l.exterior,
        core_l=core_l,
        core_r_interface=core_r.exterior,
        core_r=core_r,
        clad=clad,
        silicon=silicon,
    )

    resolutions = dict(
        core_l_interface={"resolution": 0.003, "distance": 1},
        core_r_interface={"resolution": 0.003, "distance": 1},
    )

    return mesh_from_OrderedDict(
        polygons, resolutions, filename=filename, default_resolution_max=0.3
    )


def mesh_coax(filename, radius_inner, radius_outer):
    core = shapely.Point(0, 0).buffer(radius_inner)
    isolator2 = shapely.Point(0, 0).buffer((radius_outer + radius_inner) / 4, resolution=32)
    isolator = shapely.Point(0, 0).buffer(radius_outer)

    polygons = OrderedDict(
        surface=shapely.LineString(isolator.exterior),
        core=core,
        isolator2=isolator2,
        isolator=isolator,
    )

    resolutions = dict(
        isolator2={"resolution": 0.5, "distance": 1},
        isolator={"resolution": 0.5, "distance": 1},
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=1)


if __name__ == "__main__":
    frequency = 10e9

    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh_coax(radius_inner=0.512, radius_outer=2.23039, filename="mesh.msh")
        mesh_coax(radius_inner=0.512, radius_outer=2.23039, filename=tmpdirname + "/mesh.msh")
        # mesh_waveguide_1(filename=tmpdirname + "/mesh.msh", wsim=10, hclad=1, hsi=1, wcore=1, hcore=.1, gap=.1)
        # mesh_waveguide_1(filename="mesh.msh", wsim=10, hclad=5, hsi=5, wcore=1, hcore=.1, gap=.1)
        mesh = Mesh.load(tmpdirname + "/mesh.msh")

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex)
    epsilon[basis0.get_dofs(elements="isolator")] = 1.29
    epsilon[basis0.get_dofs(elements="isolator2")] = 1.29
    epsilon[basis0.get_dofs(elements="core")] = 1 + 1j * 5.8e7 * 1e20 / frequency
    # epsilon[basis0.get_dofs(elements="silicon")] = 1.29
    # epsilon[basis0.get_dofs(elements="clad")] = 1.29
    # epsilon[basis0.get_dofs(elements="core_l")] = 1 + 1j * 5.8e7 * 1e20 / frequency
    # epsilon[basis0.get_dofs(elements="core_r")] = 1 + 1j * 5.8e7 * 1e20 / frequency
    basis0.plot(np.real(epsilon), colorbar=True).show()

    conductors = ["isolator2___isolator"]
    modes = compute_modes(
        basis0,
        epsilon,
        wavelength=scipy.constants.speed_of_light / frequency * 1e3,
        mu_r=1,
        num_modes=len(conductors),
        metallic_boundaries=True,
    )
    print("propagation constants", 1 / modes.n_effs)

    modes[0].show(modes[0].E.real, plot_vectors=True)

    from skfem import *
    from skfem.helpers import *

    @Functional(dtype=np.complex64)
    def current_form(w):
        return inner(np.array([w.n[1], -w.n[0]]), w.H)

    currents = np.zeros((len(conductors), len(modes)))

    for mode_i, mode in enumerate(modes):
        mode.show(mode.H.real, plot_vectors=True)

        (ht, ht_basis), (hz, hz_basis) = mode.basis.split(mode.H)
        for conductors_i, conductor in enumerate(conductors):
            facet_basis = ht_basis.boundary(facets=mesh.boundaries[conductor])
            current = abs(current_form.assemble(facet_basis, H=facet_basis.interpolate(ht)))
            currents[conductors_i, mode_i] = current

    characteristic_impedances = np.linalg.inv(currents).T @ np.linalg.inv(currents)
    print("characteristic impedances", characteristic_impedances)
