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
    import tempfile

    from shapely.geometry import box, Point

    from collections import OrderedDict

    electrode_left = box(9, 1, 10, 1.2)
    electrode_right = box(-9, 1, -10, 1.2)
    core = box(-10, 0, 10, 1)
    env = core.buffer(15, resolution=8)

    polygons = OrderedDict(
        electrode_left=electrode_left,
        electrode_right=electrode_right,
        core=core,
        env=env
    )

    resolutions = dict(
        core={"resolution": .1, "distance": .5}
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh_from_OrderedDict(polygons, resolutions, filename=tmpdirname + '/mesh.msh', default_resolution_max=2)
        mesh_from_OrderedDict(polygons, resolutions, filename='mesh.msh', default_resolution_max=2)
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis = Basis(mesh, ElementTriP0())

    epsilon = basis.ones()
    epsilon[basis.get_dofs(elements='core')] = 1000
    # basis.plot(epsilon).show()

    basis_u, u = solve_coulomb(basis, epsilon, {'electrode_left___core': 1, 'electrode_right___core': 0})

    basis_vec = basis_u.with_element(ElementVector(ElementTriP1()))
    u_grad = basis_vec.project(basis.interpolate(epsilon) * basis_u.interpolate(u).grad)

    x = np.linspace(np.min(basis.mesh.p[0]), np.max(basis.mesh.p[0]), 100)
    y = np.linspace(np.min(basis.mesh.p[1]), np.max(basis.mesh.p[1]), 100)

    stacked = np.dstack(np.meshgrid(x, y))
    reshaped = stacked.reshape(-1, 2)

    from shapely.prepared import prep

    env_preped = prep(env)
    inside = np.array([env_preped.contains(Point(point)) for point in reshaped])
    (u_grad_x, _), (u_grad_y, _) = basis_vec.split(u_grad)
    result = np.full_like(reshaped, np.nan)
    result[inside] = np.array(
        (basis_u.interpolator(u_grad_x)(reshaped[inside].T), basis_u.interpolator(u_grad_y)(reshaped[inside].T))).T
    result = result.reshape(stacked.shape)

    fig, ax = plt.subplots()

    for subdomain in basis.mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
        basis.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)

    basis_u.plot(u, ax=ax, shading='gouraud', colorbar=True)
    # basis_vec.plot(-u_grad, ax=ax)
    import matplotlib

    c_white = matplotlib.colors.colorConverter.to_rgba('black', alpha=0.1)
    c_black = matplotlib.colors.colorConverter.to_rgba('black', alpha=1)
    cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap', [c_white, c_black], 512)

    plt.streamplot(x, y, result[..., 0], result[..., 1], color=(np.linalg.norm(result, axis=-1)), density=5,
                   cmap=cmap_rb)
    plt.colorbar()
    plt.show()
