"""Waveguide analysis based on https://doi.org/10.1080/02726340290084012."""
import numpy as np
import scipy.linalg
import scipy.sparse.linalg

import skfem
from skfem import *
from skfem.helpers import *

mesh = skfem.Mesh.load('mesh.msh')
basis = Basis(mesh, ElementTriN0() * ElementTriP1())  # * ElementTriP1())

basis0 = basis.with_element(ElementTriP0())
epsilon = basis0.zeros()
epsilon[basis0.get_dofs(elements='Core')] = 3.4777 ** 2
epsilon[basis0.get_dofs(elements='Cladding')] = 1.444 ** 2
# basis0.plot(epsilon, colorbar=True).show()

wavelength = 1.55
k0 = 2 * np.pi / wavelength
one_over_u_r = 1


@BilinearForm
# def aform(E_t, E_z, lam, v_t, v_z, mu, w):
def aform(E_t, E_z, v_t, v_z, w):
    return one_over_u_r * curl(E_t) * curl(v_t) - k0 ** 2 * w['epsilon'] * dot(E_t, v_t) - one_over_u_r * dot(grad(E_z), v_t) + w['epsilon'] * inner(E_t, grad(v_z)) + w['epsilon'] * E_z*v_z


# @BilinearForm
# def gauge(E_t, E_z, lam, v_t, v_z, mu, w):
# set div E = 0 using a Lagrange multiplier
#    return dot(grad(lam), v_t) + dot(E_t, grad(mu))


@BilinearForm
# def bform(E_t, E_z, lam, v_t, v_z, mu, w):
def bform(E_t, E_z, v_t, v_z, w):
    return - one_over_u_r * dot(E_t, v_t)


A = aform.assemble(basis, epsilon=basis0.interpolate(epsilon))
B = bform.assemble(basis, epsilon=basis0.interpolate(epsilon))
# C = gauge.assemble(basis, epsilon=basis0.interpolate(epsilon))

# lams, xs = solve(*condense(A, B, D=basis.get_dofs()),
#                 solver=solver_eigen_scipy_sym(k=10, sigma=k0 ** 2 * 2.5 ** 2))

lams, xs = scipy.linalg.eig(A.todense(), B.todense())
# lams, xs = scipy.sparse.linalg.eigs(A, M=B, which='LR', sigma=k0 ** 2 * 2.5 ** 2, tol=.001)


idx = lams.argsort()[::-1]
#idx = idx[np.isfinite(lams[idx])]
lams = lams[idx]
xs = xs[:, idx]
xs = xs.astype(float)
print(lams.tolist())
print(np.sqrt(np.real(lams))/k0)

if __name__ == "__main__":
    print([lam for lam in np.sort(np.sqrt(np.real(lams)) / k0) if lam > 0])
    # ~2.5 is the expected effective refractive index of the mode

    idx = 0
    (Et, Etbasis), (Ez, Ezbasis), *_ = basis.split(xs[:, idx])
    print(np.sqrt(np.real(lams[idx]))/k0, np.sqrt(lams[idx])/k0)
    print(np.sum(np.abs(Et)), np.sum(np.abs(Ez)))

    Etbasis.plot(Et).show()
    Ezbasis.plot(Ez, colorbar=True).show()

    plot_basis = Etbasis.with_element(ElementVector(ElementTriP0()))
    Etp = plot_basis.project(Etbasis.interpolate(Et))
    (Etpx, plotx_basis), (Etpy, ploty_basis) = plot_basis.split(Etp)

    plotx_basis.plot(Etpx, colorbar=True, shading='gouraud').show()
    ploty_basis.plot(Etpy, colorbar=True, shading='gouraud').show()
    Ezbasis.plot(Ez, colorbar=True, shading='gouraud').show()
