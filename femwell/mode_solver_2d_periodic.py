# implementing https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-6-1053&id=312806

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
from skfem import Basis, BilinearForm, ElementDG, ElementTriP1, FacetBasis, solve
from skfem.helpers import d, grad
from skfem.io import from_meshio
from skfem.utils import mpc

from femwell.solver import solver_dense, solver_eigen_scipy_invert, solver_eigen_slepc


def solve_periodic(basis_epsilon_r, epsilon_r, k0):
    fbases = [
        FacetBasis(basis_epsilon_r.mesh, ElementTriP1(), facets="left"),
        FacetBasis(basis_epsilon_r.mesh, ElementTriP1(), facets="right"),
    ]
    assert np.all(fbases[0].default_parameters()["x"][1] == fbases[1].default_parameters()["x"][1])

    basis_vec = Basis(basis_epsilon_r.mesh, ElementTriP1() * ElementTriP1())

    @BilinearForm(dtype=np.complex64)
    def A(phi, k_phi, v, k_v, w):
        return -d(phi)[0] * d(v)[0] - d(phi)[1] * d(v)[1] + k0**2 * (w.epsilon) * phi * v

    @BilinearForm(dtype=np.complex64)
    def B(phi, k_phi, v, k_v, w):
        return 2j * grad(k_phi)[0] * v

    @BilinearForm(dtype=np.complex64)
    def C(phi, k_phi, v, k_v, w):
        return -k_phi * v

    @BilinearForm(dtype=np.complex64)
    def I_k_phi(phi, k_phi, v, k_v, w):
        return k_phi * k_v

    @BilinearForm(dtype=np.complex64)
    def I_phi(phi, k_phi, v, k_v, w):
        return phi * k_v

    A = (
        A.assemble(basis_vec, epsilon=basis_epsilon_r.interpolate(epsilon_r))
        + B.assemble(basis_vec)
        + I_k_phi.assemble(basis_vec)
    )
    f = -C.assemble(basis_vec) + I_phi.assemble(basis_vec)

    left = basis_vec.get_dofs(facets="left")
    right = basis_vec.get_dofs(facets="right")
    top = basis_vec.get_dofs(facets="top")
    bottom = basis_vec.get_dofs(facets="bottom")

    left = np.setdiff1d(left, top + bottom)
    right = np.setdiff1d(right, top + bottom)

    ks, xs = solve(
        *mpc(A, f, M=left, S=np.concatenate((right, top, bottom))),
        solver=solver_eigen_scipy_invert(
            k=20, which="LM", sigma=k0 * np.sqrt(epsilon_r.real.max())
        ),
    )
    (phis, basis_phi), (k_phis, basis_k_phi) = basis_vec.split(xs)

    return ks, basis_phi, phis


def plot_periodic(k, a, basis_phi, phi, num, ax):
    vminmax = np.max(np.abs(basis_phi.interpolate(phi)))
    for i_plot in range(num):
        phases = basis_phi.project(
            lambda x: np.exp(1j * k * (x[0] + i_plot * a)), dtype=np.complex64
        )
        phi_with_phase = basis_phi.project(
            basis_phi.interpolate(phi) * basis_phi.interpolate(phases),
            dtype=np.complex64,
        )
        mesh, z = basis_phi.refinterp(np.real(phi_with_phase), nrefs=3)
        im = ax.tripcolor(
            mesh.p[0] + i_plot * a,
            mesh.p[1],
            mesh.t.T,
            z,
            cmap="seismic",
            shading="gouraud",
            vmin=-vminmax,
            vmax=vminmax,
        )
        ax.set_aspect(1)
    return im


if __name__ == "__main__":
    from femwell.mesh import mesh_from_OrderedDict

    height = 5.76 / 2 + 5
    a = 0.010
    b = 0.78
    c = 0.2
    slab = 0.920 + 5
    pml = 3

    wavelength = 1
    k0 = 2 * np.pi / wavelength

    left = shapely.LineString([(0, y) for y in np.linspace(-height, height, 2)])
    right = shapely.LineString([(a, y) for y in np.linspace(-height, height, 2)])
    top = shapely.LineString([(x, height) for x in np.linspace(0, a, 2)])
    bottom = shapely.LineString([(x, -height) for x in np.linspace(0, a, 2)])

    box = shapely.box(0, -height, a, height)
    structure = shapely.box(0, -b / 2, a, b / 2)
    structure1 = shapely.box(0, height - slab, a, height)
    structure2 = shapely.box(0, -height + slab, a, -height)

    resolutions = {
        "structure": {"resolution": 0.1, "distance": 0.1},
        "hole": {"resolution": 0.1, "distance": 0.1},
    }

    mesh = from_meshio(
        mesh_from_OrderedDict(
            OrderedDict(
                left=left,
                right=right,
                top=top,
                bottom=bottom,
                structure=structure,
                structure1=structure1,
                structure2=structure2,
                box=box,
            ),
            resolutions=resolutions,
            filename="mesh.msh",
            default_resolution_max=0.05,
        )
    )

    basis_vec = Basis(mesh, ElementTriP1() * ElementTriP1())
    basis_epsilon_r = basis_vec.with_element(ElementDG(ElementTriP1()))

    epsilon_r = basis_epsilon_r.zeros(dtype=np.complex64) + 1.45
    epsilon_r[basis_epsilon_r.get_dofs(elements="box")] = 1.39
    epsilon_r **= 2
    basis_epsilon_r.plot(np.real(epsilon_r), ax=mesh.draw(), colorbar=True).show()
    basis_epsilon_r.plot(np.imag(epsilon_r), ax=mesh.draw(), colorbar=True).show()

    epsilon_r += basis_epsilon_r.project(
        lambda x: (0.5j) * (np.clip(np.abs(x[1]) - height + pml, 0, np.inf) / pml) ** 2,
        dtype=np.complex64,
    )

    ks, basis_phi, phis = solve_periodic(basis_epsilon_r, epsilon_r, k0)
    print(ks)

    plt.plot(np.real(ks))
    plt.plot(np.imag(ks))
    plt.show()

    for i, k in enumerate(ks):
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        plt.title(f"{k}")
        plot_periodic(k, a, basis_phi, phis[..., i], 100, ax)
        plt.show()
