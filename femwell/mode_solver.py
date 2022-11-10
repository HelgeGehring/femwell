"""Waveguide analysis based on https://doi.org/10.1080/02726340290084012."""
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg

from skfem import BilinearForm, Basis, ElementTriN0, ElementTriP0, ElementTriP1, ElementVector, Mesh, Functional
from skfem.helpers import curl, grad, dot, inner, cross


def compute_modes(basis_epsilon_r, epsilon_r, wavelength, mu_r, num_modes):
    k0 = 2 * np.pi / wavelength

    basis = basis_epsilon_r.with_element(ElementTriN0() * ElementTriP1())

    @BilinearForm(dtype=epsilon_r.dtype)
    def aform(e_t, e_z, v_t, v_z, w):
        return 1 / mu_r * curl(e_t) * curl(v_t) \
               - k0 ** 2 * w['epsilon'] * dot(e_t, v_t) \
               - 1 / mu_r * dot(grad(e_z), v_t) \
               + w['epsilon'] * inner(e_t, grad(v_z)) + w['epsilon'] * e_z * v_z

    @BilinearForm(dtype=epsilon_r.dtype)
    def bform(e_t, e_z, v_t, v_z, w):
        return - 1 / mu_r * dot(e_t, v_t)

    A = aform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))
    B = bform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))

    from petsc4py import PETSc
    from slepc4py import SLEPc

    A_ = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
    B_ = PETSc.Mat().createAIJ(size=B.shape, csr=(B.indptr, B.indices, B.data))

    eps = SLEPc.EPS().create()
    eps.setOperators(A_, B_)
    eps.getST().setType(SLEPc.ST.Type.SINVERT)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setTarget(k0 ** 2 * np.max(epsilon_r) ** 2)
    eps.setDimensions(num_modes)
    eps.solve()

    xr, wr = A_.getVecs()
    xi, wi = A_.getVecs()
    lams, xs = [], []
    for i in range(eps.getConverged()):
        lams.append(eps.getEigenpair(i, xr, xi))
        xs.append(np.array(xr) + 1j * np.array(xi))

    xs = np.array(xs, dtype=complex)
    lams = np.array(lams)
    xs[:, basis.split_indices()[1]] /= 1j * np.sqrt(lams[:, np.newaxis])  # undo the scaling E_3,new = beta * E_3

    return np.sqrt(lams) / k0, basis, xs


def calculate_hfield(basis, xs, beta):
    xs = xs.astype(complex)

    @BilinearForm(dtype=np.complex64)
    def aform(e_t, e_z, v_t, v_z, w):
        return (-1j * beta * e_t[1] + e_z.grad[1]) * v_t[0] \
               + (1j * beta * e_t[0] - e_z.grad[0]) * v_t[1] \
               + e_t.curl * v_z

    a_operator = aform.assemble(basis)

    @BilinearForm(dtype=np.complex64)
    def bform(e_t, e_z, v_t, v_z, w):
        return dot(e_t, v_t) + e_z * v_z

    b_operator = bform.assemble(basis)

    return scipy.sparse.linalg.spsolve(b_operator, a_operator @ xs) * 1j  # Don't understand the 1j yet


def calculate_overlap(basis, E_i, H_i, E_j, H_j):
    @Functional
    def overlap(w):
        return cross(np.conj(w['E_i'][0]), w['H_j'][0]) + cross(w['E_j'][0], np.conj(w['H_i'][0]))

    return overlap.assemble(basis, E_i=basis.interpolate(E_i), H_i=basis.interpolate(H_i),
                            E_j=basis.interpolate(E_j), H_j=basis.interpolate(H_j))


def plot_mode(basis, mode, plot_vectors=False, colorbar=True, title='E'):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    (et, et_basis), (ez, ez_basis) = basis.split(mode)

    if plot_vectors:
        fig, axs = plt.subplots(2, 1, subplot_kw=dict(aspect=1))
        for ax in axs:
            for subdomain in basis.mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
                basis.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
        et_basis.plot(et, ax=axs[0])
        ez_basis.plot(ez, ax=axs[1])

        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(axs[1].collections[0], cax=cax)
        return fig, axs

    plot_basis = et_basis.with_element(ElementVector(ElementTriP0()))
    et_xy = plot_basis.project(et_basis.interpolate(et))
    (et_x, et_x_basis), (et_y, et_y_basis) = plot_basis.split(et_xy)

    fig, axs = plt.subplots(3, 1, subplot_kw=dict(aspect=1))
    for ax in axs:
        for subdomain in basis.mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
            basis.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)

    for ax, component in zip(axs, 'xyz'):
        ax.set_title(f'${title}_{component}$')
    et_x_basis.plot(et_x, shading='gouraud', ax=axs[0])  # , vmin=np.min(mode), vmax=np.max(mode))
    et_y_basis.plot(et_y, shading='gouraud', ax=axs[1])  # , vmin=np.min(mode), vmax=np.max(mode))
    ez_basis.plot(ez, shading='gouraud', ax=axs[2])  # , vmin=np.min(mode), vmax=np.max(mode))

    if colorbar:
        for ax in axs:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(ax.collections[0], cax=cax)
    plt.tight_layout()

    return fig, axs


if __name__ == "__main__":
    from shapely.geometry import Polygon
    from collections import OrderedDict
    from femwell.mesh import mesh_from_polygons

    w_sim = 4 * 2
    h_clad = 1
    h_box = 1
    w_core = 0.5 * 3
    h_core = 0.22
    offset_heater = 2.2
    h_heater = .14
    w_heater = 2

    polygons = OrderedDict(
        core=Polygon([
            (-w_core / 2, -h_core / 2),
            (-w_core / 2, h_core / 2),
            (w_core / 2, h_core / 2),
            (w_core / 2, -h_core / 2),
        ]),
        clad=Polygon([
            (-w_sim / 2, -h_core / 2),
            (-w_sim / 2, -h_core / 2 + h_clad),
            (w_sim / 2, -h_core / 2 + h_clad),
            (w_sim / 2, -h_core / 2),
        ]),
        box=Polygon([
            (-w_sim / 2, -h_core / 2),
            (-w_sim / 2, -h_core / 2 - h_box),
            (w_sim / 2, -h_core / 2 - h_box),
            (w_sim / 2, -h_core / 2),
        ])
    )

    resolutions = dict(
        core={"resolution": 0.01, "distance": 1},
        heater={"resolution": 0.05, "distance": 1}
    )

    mesh_from_polygons(polygons, resolutions, filename='mesh.msh', default_resolution_max=.2)

    mesh = Mesh.load('mesh.msh')
    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros(dtype=complex)
    epsilon[basis0.get_dofs(elements='core')] = 3.4777 ** 2
    epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
    epsilon[basis0.get_dofs(elements='box')] = 1.444 ** 2
    # basis0.plot(epsilon, colorbar=True).show()

    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=6)

    print(lams)

    plot_mode(basis, np.real(xs[0]))
    plt.show()
    plot_mode(basis, np.imag(xs[0]))
    plt.show()

    xbs = calculate_hfield(basis, xs[0], -lams[0] * (2 * np.pi / 1.55))

    plot_mode(basis, np.real(xbs))
    plt.show()
    plot_mode(basis, np.imag(xbs))
    plt.show()

    integrals = np.zeros((len(lams),) * 2, dtype=complex)
    for i in range(len(lams)):
        for j in range(len(lams)):
            E_i = xs[i]
            E_j = xs[j]
            H_i = calculate_hfield(basis, E_i, -lams[i] * (2 * np.pi / 1.55))
            H_j = calculate_hfield(basis, E_j, -lams[j] * (2 * np.pi / 1.55))
            integrals[i, j] = calculate_overlap(basis, E_i, H_i, E_j, H_j)

    integrals_normalized = np.zeros((len(lams),) * 2, dtype=complex)
    for i in range(len(lams)):
        for j in range(len(lams)):
            integrals_normalized[i, j] = integrals[i, j] ** 2 / integrals[i, i] / integrals[j, j]

    plt.imshow(np.real(integrals_normalized))
    plt.colorbar()
    plt.show()
