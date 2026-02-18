# example from https://young.physics.ucsc.edu/115/quantumwell.pdf

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

from femwell.solver import solver_dense, solver_eigen_slepc

mesh = MeshLine(np.linspace(-5, 5, 201))
mesh = mesh.with_subdomains({"well": lambda p: abs(p[0]) < 0.5})

basis = Basis(mesh, ElementLineP1())
basis_potential = basis.with_element(ElementLineP0())

# Potential (eV)
potential = basis_potential.zeros()  # units seem off?
potential[basis_potential.get_dofs(elements="well")] = -8
basis_potential.plot(potential).show()


# K0 = 7.62036790878  # hbar^2/m0 in eV, and with distances in A
K0 = 1


@BilinearForm
def lhs(u, v, w):
    return 0.5 * K0 * inner(grad(u), grad(v)) + w["potential"] * inner(u, v)


@BilinearForm
def rhs(u, v, w):
    return inner(u, v)


A = lhs.assemble(basis, potential=basis_potential.interpolate(potential))
B = rhs.assemble(basis)

lams, xs = solve(*condense(A, B, D=basis.get_dofs()), solver=solver_dense(k=2, sigma=0, which="LR"))

print(lams)
basis.plot(xs[:, 0]).show()
basis.plot(xs[:, 1]).show()
