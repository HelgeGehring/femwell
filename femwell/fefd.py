import matplotlib.pyplot as plt
import numpy as np

from skfem import BilinearForm, Basis, ElementTriN1, ElementTriP0, ElementTriP1, Mesh, solve, FacetBasis
from skfem.helpers import curl, grad, dot, inner


def compute_modes(basis, basis_epsilon_r, epsilon_r, wavelength, mu_r, source):
    k0 = 2 * np.pi / wavelength
    one_over_u_r = 1 / mu_r

    @BilinearForm(dtype=np.complex64)
    def curl_form(Et, Ez, vt, vz, w):
        return one_over_u_r * inner(curl(Et), curl(vt)) + k0 ** 2 * w['epsilon'] * (inner(Et, vt)) \
               - one_over_u_r * (inner(grad(Ez), grad(vz))) + k0 ** 2 * w['epsilon'] * Ez * vz

    @BilinearForm(dtype=np.complex64)
    def div_form(Et, Ez, vt, vz, w):
        return dot(grad(Ez), vt) + dot(Et, grad(vz))

    A = curl_form.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))
    C = div_form.assemble(basis)

    return basis, solve(A + C, source)


if __name__ == '__main__':
    from collections import OrderedDict
    from shapely.geometry import Polygon, LineString
    from mesh import mesh_from_OrderedDict

    width = 4
    length = 10.5
    pml = .5

    width_wg_1 = .5
    length_wg_1 = 5
    extra_length_wg_1 = 1

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
    core_append = Polygon([
        (-width_wg_1 / 2, -length_wg_1),
        (-width_wg_1 / 2, -length_wg_1 - extra_length_wg_1),
        (width_wg_1 / 2, -length_wg_1 - extra_length_wg_1),
        (width_wg_1 / 2, -length_wg_1),
    ])

    source = LineString([
        (width_wg_1 / 2, -length_wg_1 / 2),
        (-width_wg_1 / 2, -length_wg_1 / 2)
    ])

    polygons = OrderedDict(
        source=source,
        core=core,
        core_append=core_append,
        box=core.buffer(1, resolution=4) - core,
        pml=core.buffer(2, resolution=4) - core.buffer(1, resolution=4),
    )

    resolutions = dict(
        core={"resolution": .05, "distance": 1},
        core_append={"resolution": .05, "distance": 1},
        box={"resolution": .05, "distance": 1},
    )

    mesh = mesh_from_OrderedDict(polygons, resolutions, filename='mesh.msh', default_resolution_max=.3)
    mesh = Mesh.load('mesh.msh')

    basis = Basis(mesh, ElementTriN1() * ElementTriP1())

    basis0 = basis.with_element(ElementTriP0())
    epsilon = basis0.zeros(dtype=complex) + 1.444 ** 2
    epsilon[basis0.get_dofs(elements='core')] = 2.8 ** 2
    epsilon[basis0.get_dofs(elements='core_append')] = 2.8 ** 2
    epsilon[basis0.get_dofs(elements='pml')] = (1.444 + 1j) ** 2
    basis0.plot(np.real(epsilon)).show()

    basis_source = FacetBasis(mesh, basis.elem, facets=mesh.boundaries['core___core_append'])

    source = basis_source.project(
        lambda x: [np.array([0 + 0 * x[0], 0 * np.exp(0 * x[0] ** 2) + 0 * x[0]]), 1 + 0 * x[0]])
    source = source.astype(complex)
    source *= 1j

    basis, x = compute_modes(basis, basis0, epsilon, 1.55, 1, source)

    from mode_solver import plot_mode

    plot_mode(basis, np.real(x), direction='x')
    plt.show()
    plot_mode(basis, np.imag(x), direction='x')
    plt.show()
