import scipy.linalg

from skfem import *
from skfem.helpers import *

mesh = MeshLine(np.linspace(-1, 1, 201))
mesh = mesh.with_subdomains({'core': lambda p: abs(p[0]) < .11})

basis_epsilon = Basis(mesh, ElementLineP0())

epsilon = basis_epsilon.zeros() + 1.444 ** 2
epsilon[basis_epsilon.get_dofs(elements='core')] = 3.4777 ** 2
basis_epsilon.plot(epsilon).show()

wavelength = 1.55
k0 = 2 * np.pi / wavelength


@BilinearForm
def lhs(u, v, w):
    return -1 / k0 ** 2 * dot(grad(u), grad(v)) + w['epsilon'] * inner(u, v)


@BilinearForm
def rhs(u, v, w):
    return inner(u, v)


basis = basis_epsilon.with_element(ElementLineP1())

A = lhs.assemble(basis, epsilon=basis_epsilon.interpolate(epsilon))
B = rhs.assemble(basis)

# lams, xs = scipy.linalg.eigh(A.todense(), B.todense())
lams, xs = solve(*condense(A, B, D=basis.get_dofs()), solver=solver_eigen_scipy_sym(sigma=3.55 ** 2, which='LM'))

print(np.sqrt(lams))
basis.plot(xs[:, -1]).show()
