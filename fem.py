import skfem
from skfem import Basis, ElementTriP0, ElementTriP1, ElementVector, asm, BilinearForm, solve, solver_eigen_scipy_sym, \
    condense
import numpy as np
from skfem.visuals.matplotlib import draw, plot, savefig

wavelength = 1.55e-6
k0 = 2 * np.pi / wavelength

mesh = skfem.Mesh.load('mesh.msh')
element = ElementVector(ElementTriP1(), 2)
basis = Basis(mesh, element)

basis0 = basis.with_element((ElementTriP0()))

epsilon = basis0.zeros()
epsilon[basis0.get_dofs(elements='Core')] = 3.4777 ** 2
epsilon[basis0.get_dofs(elements='Cladding')] = 1.444 ** 2


@BilinearForm
def operator(u, v, w):
    return (
            (  # top left term
                    k0 ** 2 * w['epsilon'].value * u.value[0] * v.value[0]
                    - u.grad[0][0] * v.grad[0][0]
                    - 1 / w['epsilon'].value * u.value[0] * w['epsilon'].grad[0] * v.grad[0][0]
                    - u.grad[0][1] * v.grad[0][1]
            )
            +
            (  # top right
                    - w['epsilon'].value ** (-1) * w['epsilon'].grad[1] * u.value[1] * v.grad[0][0]
            )
            +
            (  # bottom left
                    - w['epsilon'].value ** (-1) * w['epsilon'].grad[0] * u.value[0] * v.grad[1][1]
            )
            +
            (  # bottom right term
                    + k0 ** 2 * w['epsilon'].value * u.value[1] * v.value[1]
                    - u.grad[1][0] * v.grad[1][0]
                    - u.grad[1][1] * v.grad[1][1]
                    - 1 / w['epsilon'].value * u.value[1] * w['epsilon'].grad[1] * v.grad[1][1]
            )
    )


@BilinearForm
def mass(u, v, _):
    return u.value[0] * v.value[0] + u.value[1] * v.value[1]


L = asm(operator, basis, epsilon=basis0.interpolate(epsilon))
M = asm(mass, basis)

ks, u = solve(*condense(L, M, D=basis.get_dofs()), solver=solver_eigen_scipy_sym(k=8, sigma=None, which='LA'))

print(np.sqrt(ks / k0 ** 2))

plot(basis.split_bases()[0], np.real(u[basis.split_indices()[0], -1]), ax=draw(mesh), colorbar=True)
savefig('figEx.png')
plot(basis.split_bases()[1], np.real(u[basis.split_indices()[1], -1]), ax=draw(mesh), colorbar=True)
savefig('figEy.png')
plot(basis0, epsilon, ax=draw(mesh), colorbar=True)
savefig('fig.png')
