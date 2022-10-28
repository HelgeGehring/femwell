"""Waveguide analysis based on https://doi.org/10.1080/02726340290084012."""
import numpy as np

from skfem import BilinearForm, Basis, ElementTriN0, ElementTriP0, ElementTriP1, ElementVector, Mesh
from skfem.helpers import curl, grad, dot, inner


def compute_modes(basis_epsilon_r, epsilon_r, wavelength, mu_r):
    k0 = 2 * np.pi / wavelength

    basis = basis_epsilon_r.with_element(ElementTriN0() * ElementTriP1())

    @BilinearForm
    def aform(e_t, e_z, v_t, v_z, w):
        return 1 / mu_r * curl(e_t) * curl(v_t) \
               - k0 ** 2 * w['epsilon'] * dot(e_t, v_t) \
               - 1 / mu_r * dot(grad(e_z), v_t) \
               + w['epsilon'] * inner(e_t, grad(v_z)) + w['epsilon'] * e_z * v_z

    @BilinearForm
    def bform(e_t, e_z, v_t, v_z, w):
        return - 1 / mu_r * dot(e_t, v_t)

    A = aform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))
    B = bform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))

    # lams, xs = solve(*condense(A, B, D=basis.get_dofs()),
    #                solver=solver_eigen_scipy_sym(k=10, sigma=k0 ** 2 * 2.5 ** 2))
    from petsc4py import PETSc
    from slepc4py import SLEPc

    A_ = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
    B_ = PETSc.Mat().createAIJ(size=B.shape, csr=(B.indptr, B.indices, B.data))

    eps = SLEPc.EPS().create()
    eps.setOperators(A_, B_)
    eps.getST().setType(SLEPc.ST.Type.SINVERT)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setTarget(k0 ** 2 * np.max(epsilon_r) ** 2)
    eps.setDimensions(20)
    eps.solve()

    xr, xi = A_.getVecs()
    lams, xs = [], []
    for i in range(eps.getConverged()):
        lams.append(eps.getEigenpair(i, xr, xi))
        xs.append(np.array(xr) + 1j * np.array(xi))

    return np.sqrt(lams) / k0, basis, np.array(xs).T


if __name__ == "__main__":
    mesh = Mesh.load('mesh.msh')
    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros()
    epsilon[basis0.get_dofs(elements='Core')] = 3.4777 ** 2
    epsilon[basis0.get_dofs(elements='Cladding')] = 1.444 ** 2
    # basis0.plot(epsilon, colorbar=True).show()

    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1)

    print(lams)

    idx = 0
    xs = np.real(xs)
    (et, et_basis), (ez, ez_basis), *_ = basis.split(xs[:, idx])
    print(lams[idx])
    print(np.sum(np.abs(et)), np.sum(np.abs(ez)))

    et_basis.plot(et).show()
    ez_basis.plot(ez, colorbar=True).show()

    plot_basis = et_basis.with_element(ElementVector(ElementTriP0()))
    et_xy = plot_basis.project(et_basis.interpolate(et))
    (et_x, et_x_basis), (et_y, et_y_basis) = plot_basis.split(et_xy)

    et_x_basis.plot(et_x, colorbar=True, shading='gouraud').show()
    et_y_basis.plot(et_y, colorbar=True, shading='gouraud').show()
    ez_basis.plot(ez, colorbar=True, shading='gouraud').show()
