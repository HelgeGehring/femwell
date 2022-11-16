import tempfile

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np

from shapely.geometry import Polygon

import skfem
from skfem.io.meshio import from_meshio
from skfem import Mesh, Basis, ElementTriP0, ElementVector

from femwell.mode_solver import compute_modes, plot_mode
from femwell.mesh import mesh_from_OrderedDict


def mesh_waveguide(filename, wsim, hclad, hbox, wcore, hcore):
    polygons = OrderedDict(
        core=Polygon([
            (-wcore / 2, -hcore / 2),
            (-wcore / 2, hcore / 2),
            (wcore / 2, hcore / 2),
            (wcore / 2, -hcore / 2),
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
        ])
    )

    resolutions = dict(
        core={"resolution": 0.03, "distance": 1}
    )

    return mesh_from_OrderedDict(polygons, resolutions, filename=filename, default_resolution_max=.6)


if __name__ == '__main__':
    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh = mesh_waveguide(wsim=2, hclad=.7, hbox=.5, wcore=0.5, hcore=0.22, filename=tmpdirname + '/mesh.msh')
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros()
    epsilon[basis0.get_dofs(elements='core')] = 3.4777 ** 2
    epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
    epsilon[basis0.get_dofs(elements='box')] = 1.444 ** 2
    basis0.plot(epsilon, colorbar=True).show()

    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=5)

    fig, axs = plot_mode(basis, xs[0], colorbar=False)
    plt.show()
