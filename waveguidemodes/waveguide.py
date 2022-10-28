from collections import OrderedDict
import numpy as np

from shapely.geometry import Polygon
from skfem import Mesh, Basis, ElementTriP0, ElementVector

from waveguidemodes.mode_solver import compute_modes
from waveguidemodes.mesh import mesh_from_polygons
import gmsh


def mesh_waveguide(wsim, hclad, hbox, wcore, hcore):
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

    resolutions = {
        'core': {"resolution": 0.02, "distance": 1}
    }

    mesh = mesh_from_polygons(polygons, resolutions)

    gmsh.write("../mesh.msh")
    gmsh.clear()
    mesh.__exit__()


mesh_waveguide(wsim=2, hclad=1, hbox=1, wcore=0.5, hcore=0.22)

mesh = Mesh.load('../mesh.msh')
basis0 = Basis(mesh, ElementTriP0(), intorder=4)
epsilon = basis0.zeros()
epsilon[basis0.get_dofs(elements='core')] = 3.4777 ** 2
epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
epsilon[basis0.get_dofs(elements='box')] = 1.444 ** 2
basis0.plot(epsilon, colorbar=True).show()

lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1)

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
