import matplotlib.pyplot as plt
import numpy as np

from skfem import *
from skfem.helpers import *

from mesh import mesh_from_OrderedDict


def solve_coulomb(basis_epsilon_r, epsilon_r, fixed_boundaries):
    basis = Basis(basis_epsilon_r.mesh, ElementTriP1())

    @BilinearForm
    def coulomb(u, v, w):
        return w['epsilon_r'] * inner(grad(u), grad(v))

    A = coulomb.assemble(basis,
                         epsilon_r=basis_epsilon_r.interpolate(epsilon_r))

    u = basis.zeros()
    for key, value in fixed_boundaries.items():
        u[basis.get_dofs(key)] = value

    return basis, solve(*condense(A, x=u, D={key: basis.get_dofs(key) for key in fixed_boundaries}))


if __name__ == '__main__':
    # Reproduce https://doi.org/10.3390/

    import tempfile

    from shapely.geometry import box

    from collections import OrderedDict

    electrode_start_x = 1.532 / 2 + 2.629
    electrode_width = 4.4

    electrode_left = box(-electrode_start_x - electrode_width, .5, - electrode_start_x, .5 + 1.8)
    electrode_right = box(electrode_start_x, .5, electrode_start_x + electrode_width, .5 + 1.8)
    slab = box(-10, 0, 10, .5)
    core = box(-1.532 / 2, .5, 1.532 / 2, .8)
    env = slab.buffer(20, resolution=8)

    polygons = OrderedDict(
        electrode_left=electrode_left,
        electrode_right=electrode_right,
        slab=slab,
        core=core,
        env=env
    )

    resolutions = dict(
        slab={"resolution": .02, "distance": .1},
        core={"resolution": .02, "distance": .1}
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh_from_OrderedDict(polygons, resolutions, filename=tmpdirname + '/mesh.msh', default_resolution_max=5)
        mesh_from_OrderedDict(polygons, resolutions, filename='mesh.msh', default_resolution_max=5)
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis = Basis(mesh, ElementTriP1())
    basis_epsilon_r = basis.with_element(ElementTriP0())

    epsilon = basis_epsilon_r.zeros() + 3.9
    epsilon[basis_epsilon_r.get_dofs(elements='slab')] = 128.4
    epsilon[basis_epsilon_r.get_dofs(elements='core')] = 7.5
    # basis.plot(epsilon).show()

    basis_u, u = solve_coulomb(basis_epsilon_r, epsilon, {'electrode_left___slab': 1, 'electrode_right___slab': 0})

    fig, ax = plt.subplots()

    for subdomain in basis_epsilon_r.mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
        basis_epsilon_r.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)

    basis_u.plot(u, ax=ax, shading='gouraud', colorbar=True)
    # basis_vec.plot(-u_grad, ax=ax)
    plt.show()

    fig, ax = plt.subplots()

    for subdomain in basis_epsilon_r.mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
        basis_epsilon_r.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
    basis_grad = basis_u.with_element(ElementDG(basis_u.elem))
    basis_u.plot(basis_u.project(basis_epsilon_r.interpolate(epsilon) * basis_u.interpolate(u).grad[0]),
                 ax=ax, shading='gouraud', colorbar=True)
    plt.show()
