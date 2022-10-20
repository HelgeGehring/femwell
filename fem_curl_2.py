import matplotlib.pyplot as plt

import skfem
from skfem import Basis, ElementTriP0, ElementTriP1, ElementTriP2, ElementTriP3, ElementTriN0, ElementVector, asm, \
    BilinearForm, \
    solve, \
    solver_eigen_scipy_sym, solver_eigen_scipy, \
    condense, ElementComposite
from skfem.helpers import curl, dot, ddot, grad, inner
import numpy as np
from skfem.visuals.matplotlib import draw, plot, savefig

# wavelength = 1.55e-6
# k0 = 2 * np.pi / wavelength

mesh = skfem.MeshTri().init_tensor(np.linspace(0, 1, 30), np.linspace(0, .5, 30))
basis = Basis(mesh, ElementComposite(ElementTriN0(), ElementTriP1()))

basis0 = basis.with_element(ElementTriP0())

epsilon = basis0.zeros() + 1
one_over_u_r = 1

plot(basis0, epsilon, ax=draw(mesh), colorbar=True)
plt.show()


@BilinearForm
def s_tt_ij(N_i, L_i, N_j, L_j, w):
    return one_over_u_r * inner(curl(N_i), curl(N_j))


@BilinearForm
def t_tt_ij(N_i, L_i, N_j, L_j, w):
    return w['epsilon'] * inner(N_i, N_j)


@BilinearForm
def s_zz_ij(N_i, L_i, N_j, L_j, w): \
        return one_over_u_r * inner(grad(L_i), grad(L_j))


@BilinearForm
def t_zz_ij(N_i, L_i, N_j, L_j, w):
    return w['epsilon'] * L_i * L_j


s_ij = asm(s_tt_ij, basis, epsilon=basis0.interpolate(epsilon)) \
       + asm(s_zz_ij, basis, epsilon=basis0.interpolate(epsilon))
t_ij = asm(t_tt_ij, basis, epsilon=basis0.interpolate(epsilon)) \
       + asm(t_zz_ij, basis, epsilon=basis0.interpolate(epsilon))

eigenValues, eigenVectors = solve(*condense(s_ij, t_ij, D=mesh.boundary_nodes()),
                                  solver=solver_eigen_scipy(k=6, sigma=1, which='LM'))

idx = eigenValues.argsort()[::-1]
eigenValues = eigenValues[idx]
eigenVectors = eigenVectors[:, idx]

print(eigenValues)

plot(basis.split_bases()[0], eigenVectors[basis.split_indices()[0], 0], ax=draw(mesh, boundaries_only=True),
     colorbar=True)
plt.show()

plot(basis.split_bases()[1], eigenVectors[basis.split_indices()[1], 0], ax=draw(mesh, boundaries_only=True),
     colorbar=True)
plt.show()

print(np.max(eigenVectors[basis.split_indices()[0]]), np.max(eigenVectors[basis.split_indices()[1]]))
