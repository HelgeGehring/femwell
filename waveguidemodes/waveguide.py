import tempfile

from collections import OrderedDict
import numpy as np

from shapely.geometry import Polygon

import skfem
from skfem.io.meshio import from_meshio
from skfem import Mesh, Basis, ElementTriP0, ElementVector

from waveguidemodes.mode_solver import compute_modes
from waveguidemodes.mesh import mesh_from_polygons


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
        core={"resolution": 0.05, "distance": 10},
        clad={"resolution": 0.1, "distance": 10},
        box={"resolution": 0.1, "distance": 10}
    )

    return mesh_from_polygons(polygons, resolutions, filename=filename, default_resolution_max=.1)


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

    idx = 0
    xs = np.real(xs)
    (et, et_basis), (ez, ez_basis), *_ = basis.split(xs[:, idx])
    plot_basis = et_basis.with_element(ElementVector(ElementTriP0()))
    et_xy = plot_basis.project(et_basis.interpolate(et))
    (et_x, et_x_basis), (et_y, et_y_basis) = plot_basis.split(et_xy)

    print(lams)
    print(np.sum(np.abs(et_x)), np.sum(np.abs(et_y)), np.sum(np.abs(ez)))

    et_x_basis.plot(et_x, colorbar=True, shading='gouraud').show()
    et_y_basis.plot(et_y, colorbar=True, shading='gouraud').show()
    ez_basis.plot(ez, colorbar=True, shading='gouraud').show()

    print(np.sum(np.abs(xs[basis.get_dofs(), idx])) / np.sum(np.abs(xs[:, idx])))
