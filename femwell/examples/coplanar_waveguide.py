import tempfile

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import shapely.ops

from shapely.geometry import LineString, box
from skfem import Mesh, Basis, ElementTriP0

from femwell.mode_solver import compute_modes, plot_mode, calculate_hfield
from femwell.mesh import mesh_from_OrderedDict


def mesh_waveguide(filename, wsim, hclad, hsi, wcore, hcore):
    core = box(-wcore / 2, -hcore / 2, wcore / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)

    all = shapely.ops.unary_union([core, clad, silicon])

    polygons = OrderedDict(
        boundary=LineString(all.exterior),
        interface=LineString([
            (-wsim / 2, -hcore / 2),
            (wsim / 2, -hcore / 2),
        ]),
        core_interface=core.exterior,
        core=core,
        clad=clad,
        silicon=silicon,
    )

    resolutions = dict(
        core_interface={"resolution": .003, "distance": 1},
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=1)


def mesh_waveguide_1(filename, wsim, hclad, hsi, wcore, hcore, gap):
    core_l = box(-wcore - gap, -hcore / 2, -gap / 2, hcore / 2)
    core_r = box(wcore + gap, -hcore / 2, gap / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)

    all = shapely.ops.unary_union([core_l, core_r, clad, silicon])

    polygons = OrderedDict(
        boundary=LineString(all.exterior),
        interface=LineString([
            (-wsim / 2, -hcore / 2),
            (wsim / 2, -hcore / 2),
        ]),
        core_l_interface=core_l.exterior,
        core_l=core_l,
        core_r_interface=core_r.exterior,
        core_r=core_r,
        clad=clad,
        silicon=silicon,
    )

    resolutions = dict(
        core_l_interface={"resolution": .003, "distance": 1},
        core_r_interface={"resolution": .003, "distance": 1},
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=1)


if __name__ == '__main__':
    omega = 9e9
    print('lambda ', 2 * np.pi * scipy.constants.c * 1e3 / omega)

    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh_waveguide_1(wsim=10, hclad=4, hsi=1, wcore=2, hcore=.02, gap=1,
                       filename=tmpdirname + '/mesh.msh')
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex)
    epsilon[basis0.get_dofs(elements='clad')] = 1
    # epsilon[basis0.get_dofs(elements='core')] = + 1j * 5.8e7 * 1e20 / omega
    epsilon[basis0.get_dofs(elements='core_l')] = + 1j * 5.8e7 * 1e20 / omega
    epsilon[basis0.get_dofs(elements='core_r')] = + 1j * 5.8e7 * 1e20 / omega
    epsilon[basis0.get_dofs(elements='silicon')] = 4
    basis0.plot(np.real(epsilon), colorbar=True).show()

    print(2 * np.pi * scipy.constants.c * 1e3 / omega, '<--')
    print(2 * np.pi / omega)
    conductors = ['core_l_interface', 'core_r_interface']
    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1/.1, mu_r=1,
                                    num_modes=len(conductors), metallic_boundaries=True)
    print('lams', lams)

    fig, axs = plot_mode(basis, np.real(xs[1]), plot_vectors=True)
    plt.show()

    # plot_mode(basis, np.real(xbs), plot_vectors=True)
    # plt.show()

    from skfem import *
    from skfem.helpers import *


    @Functional
    def current_form(w):
        return inner(np.array([w.n[1], -w.n[0]]), w.H) * 1e-6


    currents = np.zeros((len(conductors), len(lams)))

    for mode_i in range(len(lams)):
        xbs = calculate_hfield(basis, xs[mode_i], lams[mode_i] / omega)
        (ht, ht_basis), (hz, hz_basis) = basis.split(xbs)
        for conductors_i, conductor in enumerate(conductors):
            facet_basis = FacetBasis(ht_basis.mesh, ht_basis.elem, facets=mesh.boundaries[conductor])
            current = abs(current_form.assemble(facet_basis, H=facet_basis.interpolate(ht)))
            print(f'mode {mode_i} current in ' + conductor + '\t', current)
            currents[conductors_i, mode_i] = current

    print('currents', currents)
    print(np.linalg.inv(currents))

    characteristic_impedances = np.linalg.inv(currents).T @ np.linalg.inv(currents)
    print(characteristic_impedances)
