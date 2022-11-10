import matplotlib.pyplot as plt

from skfem import *
from skfem.helpers import *


def solve_coulomb(basis_epsilon_r, epsilon_r, fixed_boundaries):
    basis = basis_epsilon_r.with_element(ElementTriP1())

    @BilinearForm
    def coulomb(u, v, w):
        return w['epsilon_r'] * inner(div(u), div(v))

    A = coulomb.assemble(basis,
                         epsilon_r=basis_epsilon_r.interpolate(epsilon_r))

    u = basis.zeros()
    for key, value in fixed_boundaries.items():
        u[basis.get_dofs(key)] = value

    return basis, solve(*condense(A, x=u, D={key: basis.get_dofs(key) for key in fixed_boundaries}))


if __name__ == '__main__':
    mesh = MeshTri().refined(4)
    basis = Basis(mesh, ElementTriP1())

    epsilon = basis.project(lambda x: 1 + 5 * ((x[0] - .5) ** 2 + (x[1] - .5) ** 2 < .2 ** 2))

    basis_u, u = solve_coulomb(basis, epsilon, {'left': 1, 'right': 0})

    basis_vec = basis_u.with_element(ElementVector(ElementTriP1()))
    u_grad = basis_vec.project(basis_u.interpolate(u).grad)
    fig, ax = plt.subplots()
    # ax = basis.draw()
    basis_u.plot(u, ax=ax, shading='gouraud', colorbar=True)
    basis_vec.plot(-u_grad, ax=ax)
    plt.show()
