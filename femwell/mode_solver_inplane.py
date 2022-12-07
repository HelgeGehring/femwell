import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg

from skfem import BilinearForm, Basis, ElementTriN1, ElementTriP0, ElementTriP1, ElementVector, Mesh, Functional, solve, \
    condense, solver_eigen_scipy_sym, solver_eigen_scipy, FacetBasis, asm, ElementTriN2, ElementTriP2
from skfem.helpers import curl, grad, dot, inner, cross
from femwell.mode_solver import solver_slepc


def compute_modes(basis_epsilon_r, epsilon_r, wavelength, mu_r, num_modes, phase_x):
    k0 = 2 * np.pi / wavelength
    one_over_u_r = 1

    basis = basis_epsilon_r.with_element(ElementTriN2() * ElementTriP2())

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

    @BilinearForm(dtype=complex)
    def penalty(u, u1, v, v1, w):
        u1 = (w.idx[0] == 0) * u
        u2 = (w.idx[0] == 1) * u
        v1 = (w.idx[1] == 0) * v
        v2 = (w.idx[1] == 1) * v
        ju = u1 - w['phase'] * u2
        jv = v1 - w['phase'] * v2

        u1 = (w.idx[0] == 0) * u1
        u2 = (w.idx[0] == 1) * u1
        v1 = (w.idx[1] == 0) * v1
        v2 = (w.idx[1] == 1) * v1
        ju_ = u1 - w['phase'] * u2
        jv_ = v1 - w['phase'] * v2

        return 1. / 1e-2 * inner(ju, jv) + 1. / 1e-2 * inner(ju_, jv_)

    fbases = [
        FacetBasis(basis.mesh, basis.elem, facets='left'),
        FacetBasis(basis.mesh, basis.elem, facets='right'),
    ]
    D1 = asm(penalty, fbases, fbases, phase=phase_x)

    fbases = [
        FacetBasis(basis.mesh, basis.elem, facets='top'),
        FacetBasis(basis.mesh, basis.elem, facets='bottom'),
    ]
    D2 = asm(penalty, fbases, fbases, phase=1)

    lams, xs = solve(*condense(A + C + D1 + D2, B, D=basis.get_dofs(), x=basis.zeros(dtype=complex)),
                     solver=solver_eigen_scipy(k=num_modes, which='LM', ))

    return np.sqrt(lams) / k0, basis, xs


if __name__ == '__main__':
    from collections import OrderedDict
    from shapely.geometry import Polygon, box, Point, LineString
    from mesh import mesh_from_OrderedDict

    width = 1
    cell_width = .65
    cell_height = 5

    structure = box(0, -width / 2, cell_width, width / 2)
    hole = structure.centroid.buffer(.2)
    cell = box(0, -cell_height / 2, cell_width, cell_height / 2) - structure

    polygons = OrderedDict(
        left=LineString(((0, -cell_height / 2), (0, cell_height / 2))),
        right=LineString(((cell_width, -cell_height / 2), (cell_width, cell_height / 2))),
        top=LineString(((0, cell_height / 2), (cell_width, cell_height / 2))),
        bottom=LineString(((0, -cell_height / 2), (cell_width, -cell_height / 2))),
        # hole=hole,
        structure=structure,
        cell=cell.geoms[0],
        cell1=cell.geoms[1]
    )

    resolutions = dict(
        core={"resolution": .02, "distance": 1},
    )

    mesh = mesh_from_OrderedDict(polygons, resolutions, filename='mesh.msh', default_resolution_max=.05)
    mesh = Mesh.load('mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros(dtype=complex) + 1.444 ** 2
    epsilon[basis0.get_dofs(elements='structure')] = 2.8 ** 2

    phases = np.linspace(0, 1, 20)
    results = []
    for phase in phases:
        lams, basis, xs = compute_modes(basis0, epsilon, 1.55, 1, 5, phase_x=np.exp(1j * np.pi * phase))
        print(lams)
        results.append((lams, basis, xs))

    lams = np.array([result[0] for result in results])
    plt.plot(phases, lams)
    plt.show()

    from mode_solver import plot_mode

    plot_mode(basis, np.real(xs[:, -1]), direction='x', colorbar=True)
    plt.show()

    plot_mode(basis, np.real(xs[:, -1]), direction='x', colorbar='same')
    plt.show()
