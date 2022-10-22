"""Waveguide cutoff analysis."""

import numpy as np
import skfem
from skfem import *
from skfem.helpers import *

mesh = skfem.Mesh.load('mesh.msh')
basis = Basis(mesh, ElementTriN0() * ElementTriP1())  # * ElementTriP1())

basis0 = basis.with_element(ElementTriP0())
epsilon = basis0.zeros()
epsilon[basis0.get_dofs(elements='Core')] = 3.4777 ** 2
epsilon[basis0.get_dofs(elements='Cladding')] = 1.444 ** 2
basis0.plot(epsilon, colorbar=True).show()

wavelength = 1.55
k0 = 2 * np.pi / wavelength
one_over_u_r = 1


@BilinearForm
# def aform(E_t, E_z, lam, v_t, v_z, mu, w):
def aform(E_t, E_z, v_t, v_z, w):
    return one_over_u_r * curl(E_t) * curl(v_t) - k0 ** 2 * w['epsilon'] * dot(E_t, v_t)


@BilinearForm
def gauge(E_t, E_z, lam, v_t, v_z, mu, w):
    # set div E = 0 using a Lagrange multiplier
    return dot(grad(lam), v_t) + dot(E_t, grad(mu)) + lam * v_z + E_z * mu


@BilinearForm
# def bform(E_t, E_z, lam, v_t, v_z, mu, w):
def bform(E_t, E_z, v_t, v_z, w):
    return - (one_over_u_r * dot(grad(E_z) + E_t, grad(v_z) + v_t) - k0 ** 2 * w['epsilon'] * E_z * v_z)


A = aform.assemble(basis, epsilon=basis0.interpolate(epsilon))
B = bform.assemble(basis, epsilon=basis0.interpolate(epsilon))
# C = gauge.assemble(basis, epsilon=basis0.interpolate(epsilon))

lams, xs = solve(*condense(A, B, D=basis.get_dofs()),
                 solver=solver_eigen_scipy_sym(k=5, sigma=k0 ** 2 * 2.5 ** 2))

if __name__ == "__main__":
    print(np.sqrt(lams) / k0)

    (Et, Etbasis), (Ez, Ezbasis), *_ = basis.split(xs[:, 0])
    print(np.max(np.abs(Et)), np.max(np.abs(Ez)))

    Etbasis.plot(Et).show()
    Ezbasis.plot(Ez, colorbar=True).show()
