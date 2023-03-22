import numpy as np
from skfem import (
    Basis,
    BilinearForm,
    ElementLineP0,
    ElementLineP1,
    MeshLine,
    condense,
    solve,
)
from skfem.helpers import dot, grad, inner
from skfem.utils import solver_eigen_scipy

from femwell.solver import solver_eigen_slepc

mesh = MeshLine(np.linspace(-1, 1, 201))
mesh = mesh.with_subdomains({"well": lambda p: abs(p[0]) < 0.2})

basis = Basis(mesh, ElementLineP1())
basis_potential = basis.with_element(ElementLineP0())

# Potential (eV)
potential = basis_potential.zeros() + 1000  # units seem off?
potential[basis_potential.get_dofs(elements="well")] = 0
basis_potential.plot(potential).show()


K0 = 7.62036790878  # hbar^2/m0 in eV, and with distances in A


@BilinearForm
def lhs(u, v, w):
    return K0 * dot(grad(u), grad(v)) + w["potential"] * inner(u, v)


@BilinearForm
def rhs(u, v, w):
    return inner(u, v)


A = lhs.assemble(basis, potential=basis_potential.interpolate(potential))
B = rhs.assemble(basis)

lams, xs = solve(
    *condense(A, B, D=basis.get_dofs()), solver=solver_eigen_scipy(k=2, sigma=0, which="LM")
)

print(np.sqrt(lams))
basis.plot(xs[:, 0]).show()
basis.plot(xs[:, 1]).show()
