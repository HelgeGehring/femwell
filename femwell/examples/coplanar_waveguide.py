import tempfile

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely.ops

from shapely.geometry import LineString, box
from skfem import Mesh, Basis, ElementTriP0

from femwell.mode_solver import compute_modes, plot_mode, calculate_hfield, calculate_overlap
from femwell.mesh import mesh_from_OrderedDict


def mesh_waveguide(filename, wsim, hclad, hsi, wcore, hcore, gap):
    core = box(-wcore / 2, -hcore / 2, wcore / 2, hcore / 2)
    core_l = box(-wcore / 2 - gap, -hcore / 2, -wsim / 2, hcore / 2)
    core_r = box(wcore / 2 + gap, -hcore / 2, wsim / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)

    all = shapely.ops.unary_union([core, core_l, core_r, clad, silicon])

    polygons = OrderedDict(
        boundary=LineString(all.exterior),
        interface=LineString([
            (-wsim / 2, -hcore / 2),
            (wsim / 2, -hcore / 2),
        ]),
        core_interface=core.exterior,
        core=core,
        core_l_interface=core_l.exterior,
        core_l=core_l,
        core_r_interface=core_r.exterior,
        core_r=core_r,
        clad=clad,
        silicon=silicon,
    )

    resolutions = dict(
        interface={"resolution": 1, "distance": 10},
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=10)


if __name__ == '__main__':
    omega = 1e-5

    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh_waveguide(wsim=200, hclad=50, hsi=50, wcore=10, hcore=1, gap=20,
                       filename=tmpdirname + '/mesh.msh')
        mesh_waveguide(wsim=200, hclad=50, hsi=50, wcore=10, hcore=1, gap=20,
                       filename='mesh.msh')
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex)
    epsilon[basis0.get_dofs(elements='core')] = 1 + 1j * 5.8e7 / omega
    epsilon[basis0.get_dofs(elements='clad')] = 1
    epsilon[basis0.get_dofs(elements='core_r')] = 1 + 1j * 5.8e7 / omega
    epsilon[basis0.get_dofs(elements='core_l')] = 1 + 1j * 5.8e7 / omega
    epsilon[basis0.get_dofs(elements='silicon')] = 11.7
    # basis0.plot(np.real(epsilon), colorbar=True).show()

    print(2 * np.pi / omega)
    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=2 * np.pi / omega, mu_r=1, num_modes=5)

    # fig, axs = plot_mode(basis, np.real(xs[1]), plot_vectors=True)
    # plt.show()

    xbs = calculate_hfield(basis, xs[1], lams[1] / omega)

    # plot_mode(basis, np.real(xbs), plot_vectors=True)
    # plt.show()

    (ht, ht_basis), (hz, hz_basis) = basis.split(xbs)

    from skfem import *
    from skfem.helpers import *

    facet_basis = FacetBasis(ht_basis.mesh, ht_basis.elem, facets=mesh.boundaries['core_interface'])
    facet_basis_l = FacetBasis(ht_basis.mesh, ht_basis.elem, facets=mesh.boundaries['core_l_interface'])
    facet_basis_r = FacetBasis(ht_basis.mesh, ht_basis.elem, facets=mesh.boundaries['core_r_interface'])
    facet_basis_b = FacetBasis(ht_basis.mesh, ht_basis.elem, facets=mesh.boundaries['boundary'])


    @Functional
    def current(w):
        return inner(np.array([w.n[1], -w.n[0]]), w.H)


    print('current', abs(current.assemble(facet_basis, H=facet_basis.interpolate(ht))))
    print('current', abs(current.assemble(facet_basis_l, H=facet_basis_l.interpolate(ht))))
    print('current', abs(current.assemble(facet_basis_r, H=facet_basis_r.interpolate(ht))))
    print('current', abs(current.assemble(facet_basis_b, H=facet_basis_b.interpolate(ht))))
