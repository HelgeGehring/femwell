import tempfile

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np

from shapely.geometry import Polygon

import skfem
from skfem.io.meshio import from_meshio
from skfem import Mesh, Basis, ElementTriP0, ElementVector

from waveguidemodes.mode_solver import compute_modes, plot_mode
from waveguidemodes.mesh import mesh_from_polygons


def mesh_waveguide(filename, wsim, hclad, hbox, wcore, hcore, core_offset):
    polygons = OrderedDict(
        core=Polygon([
            (-wcore / 2 - core_offset / 2, -hcore / 2),
            (-wcore / 2 - core_offset / 2, hcore / 2),
            (wcore / 2 - core_offset / 2, hcore / 2),
            (wcore / 2 - core_offset / 2, -hcore / 2),
        ]),
        core2=Polygon([
            (-wcore / 2 + core_offset / 2, -hcore / 2),
            (-wcore / 2 + core_offset / 2, hcore / 2),
            (wcore / 2 + core_offset / 2, hcore / 2),
            (wcore / 2 + core_offset / 2, -hcore / 2),
        ]),
        clad=Polygon([
            (-wsim / 2, -hcore / 2),
            (-wsim / 2, -hcore / 2 + hclad),
            (wsim / 2, -hcore / 2 + hclad),
            (wsim / 2, -hcore / 2),
        ]),
        box=Polygon([
            (-wsim / 2, -hcore / 2),
            (-wsim / 2, -hcore / 2 - hbox),
            (wsim / 2, -hcore / 2 - hbox),
            (wsim / 2, -hcore / 2),
        ]),
    )

    resolutions = dict(
        core={"resolution": .3, "distance": 10},
        core2={"resolution": .3, "distance": 10},
        clad={"resolution": .3, "distance": 10},
        box={"resolution": .3, "distance": 10}
    )

    return mesh_from_polygons(polygons, resolutions, filename=filename, default_resolution_max=1)


if __name__ == '__main__':
    omega = 1e-5

    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh = mesh_waveguide(wsim=50, hclad=26, hbox=26, wcore=10, hcore=1, core_offset=20,
                              filename=tmpdirname + '/mesh.msh')
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex)
    epsilon[basis0.get_dofs(elements='core')] = 1 + 1j * 5.8e7 / omega
    epsilon[basis0.get_dofs(elements='clad')] = 1
    epsilon[basis0.get_dofs(elements='core2')] = 1 + 1j * 5.8e7 / omega
    epsilon[basis0.get_dofs(elements='box')] = 11.7
    basis0.plot(np.real(epsilon), colorbar=True).show()

    print(2 * np.pi / omega)
    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=2 * np.pi / omega, mu_r=1, num_modes=5)

    fig, axs = plot_mode(basis, xs[0], plot_vectors=True)
    plt.show()
