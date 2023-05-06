from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
from skfem import *
from skfem.helpers import grad, inner
from skfem.io import from_meshio

from femwell.mesh import mesh_from_OrderedDict

waveguide = shapely.box(0, 0, 10, 1)
mesh = from_meshio(
    mesh_from_OrderedDict(
        OrderedDict(waveguide=waveguide), resolutions={}, default_resolution_max=0.1
    )
)

basis = Basis(mesh, ElementTriP1())
basis0 = basis.with_element(ElementTriP0())
epsilon_r = basis0.zeros(dtype=complex) + 1
dofs = basis0.get_dofs(elements=lambda x: x[0] > 5)
epsilon_r[dofs] += basis0.project(lambda x: np.maximum(0, x[0] - 5) ** 2 * 0.05j, dtype=complex)[
    dofs
]
basis0.plot(epsilon_r.imag, shading="gouraud", colorbar=True)
plt.show()

mu_r = 1
k0 = 4


@BilinearForm(dtype=complex)
def maxwell(u, v, w):
    return 0.5 * (1 / w.epsilon_r * inner(grad(u), grad(v))) - k0**2 * inner(u, v)


A = maxwell.assemble(basis, epsilon_r=basis0.interpolate(epsilon_r))
B = basis.zeros(dtype=complex)
dofs = basis.get_dofs(lambda x: x[0] == np.min(x[0]))
B[dofs] = basis.project(lambda x: np.sin(np.pi * x[1]), dtype=complex)[dofs]

C = solve(
    *condense(
        A,
        x=B,
        D=basis.get_dofs(
            {
                lambda x: x[0] == np.min(x[0]),
                lambda x: x[1] == np.min(x[1]),
                lambda x: x[1] == np.max(x[1]),
            }
        ),
        expand=True,
    )
)
# C = solve(A,B)

basis.plot(C.real, shading="gouraud", colorbar=True)
plt.show()
print(C)
