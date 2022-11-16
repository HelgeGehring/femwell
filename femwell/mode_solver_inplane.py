import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg

from skfem import BilinearForm, Basis, ElementTriN1, ElementTriP0, ElementTriP1, ElementVector, Mesh, Functional, solve, \
    condense, solver_eigen_scipy_sym, solver_eigen_scipy
from skfem.helpers import curl, grad, dot, inner, cross


def compute_modes(basis_epsilon_r, epsilon_r, wavelength, mu_r, num_modes):
    k0 = 2 * np.pi / wavelength
    one_over_u_r = 1

    basis = basis_epsilon_r.with_element(ElementTriN1() * ElementTriP1())

    @BilinearForm(dtype=complex)
    def aform(E, lam, v, mu, w):
        return one_over_u_r * curl(E) * curl(v)

    @BilinearForm(dtype=complex)
    def gauge(E, lam, v, mu, w):
        # set div E = 0 using a Lagrange multiplier
        return dot(grad(lam), v) + dot(E, grad(mu))

    @BilinearForm(dtype=complex)
    def bform(E, lam, v, mu, w):
        return w['epsilon'] * dot(E, v)

    A = aform.assemble(basis)
    B = bform.assemble(basis, epsilon=basis0.interpolate(epsilon_r))
    C = gauge.assemble(basis)

    lams, xs = solve(*condense(A + C, B, D=basis.get_dofs(), x=basis.zeros(dtype=complex)),
                     solver=solver_eigen_scipy_sym(k=1))

    return np.sqrt(lams) / k0, basis, xs


if __name__ == '__main__':
    from collections import OrderedDict
    from shapely.geometry import Polygon
    from mesh import mesh_from_OrderedDict

    width = 4
    length = 10.5
    pml = .5

    width_wg_1 = .5
    length_wg_1 = 5

    width_wg_2 = 2
    length_wg_2 = 5

    core = Polygon([
        (-width_wg_1 / 2, -length_wg_1),
        (-width_wg_1 / 2, 0),
        (-width_wg_2 / 2, 0),
        (-width_wg_2 / 2, length_wg_2),
        (width_wg_2 / 2, length_wg_2),
        (width_wg_2 / 2, 0),
        (width_wg_1 / 2, 0),
        (width_wg_1 / 2, -length_wg_1),
    ])

    polygons = OrderedDict(
        core=core,
        box=core.buffer(1, resolution=4)-core,
        pml=core.buffer(2, resolution=4)-core.buffer(1, resolution=4),
    )

    resolutions = dict(
        core={"resolution": .02, "distance": 1},
    )

    mesh = mesh_from_OrderedDict(polygons, resolutions, filename='mesh.msh', default_resolution_max=.3)
    mesh = Mesh.load('mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros(dtype=complex) + 1.444 ** 2
    epsilon[basis0.get_dofs(elements='core')] = 2.8 ** 2
    epsilon[basis0.get_dofs(elements='pml')] = (1.444 + 1j) ** 2

    lams, basis, xs = compute_modes(basis0, epsilon, 1.55, 1, 1)

    from mode_solver import plot_mode

    plot_mode(basis, np.real(xs[:, 0]), direction='y')
    plt.show()
