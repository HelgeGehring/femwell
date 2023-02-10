# implementing https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-6-1053&id=312806

import numpy as np
from skfem import Basis, BilinearForm, ElementTriP1, FacetBasis, solve
from skfem.helpers import d, grad

from femwell.solver import solver_dense, solver_eigen_scipy_invert, solver_eigen_slepc
from femwell.utils import mpc_symmetric


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
        *mpc_symmetric(A, f, M=left, S=np.concatenate((right, top, bottom))),
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
