# formula from https://doi.org/10.1114/1.1352640

from collections import OrderedDict

import matplotlib
import numpy as np
from shapely import LineString, box
from skfem import (
    Basis,
    BilinearForm,
    ElementDG,
    ElementTriP0,
    ElementTriP3,
    MeshTri,
    condense,
    solve,
)
from skfem.helpers import grad, inner
from skfem.io import from_meshio

from femwell.mesh import mesh_from_OrderedDict

mesh = MeshTri().refined(3).scaled([4, 1])

core = box(-5, 0, 5, 10) - box(-2, 0, 2, 8).buffer(1)

polygons = OrderedDict(
    left=LineString([(-5, 0), (-2, 0)]),
    right=LineString([(2, 0), (5, 0)]),
    core=core,
    box=core.buffer(10),
)

resolutions = dict(
    core={"resolution": 0.3, "distance": 1},
)

mesh = from_meshio(
    mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=2, filename="mesh.msh")
)

mesh.draw().show()

basis = Basis(mesh, ElementTriP3())
basis_current = Basis(mesh, ElementDG(ElementTriP3()))
basis0 = basis.with_element(ElementTriP0())


@BilinearForm
def laplace_equation(u, v, w):
    return -inner(grad(u), w["sigma"] * grad(v))


sigma = basis0.project(lambda x: 5 * (x[0] < 2) + 1)
sigma = basis0.zeros() + 0.1
sigma[basis0.get_dofs(elements="core")] = 10
basis0.plot(sigma).show()

A = laplace_equation.assemble(basis, sigma=basis0.interpolate(sigma))

left = basis.get_dofs("left")
right = basis.get_dofs("right")
voltages = basis.zeros()
voltages[left] = 1

result = solve(*condense(A, basis.zeros(), D=left + right, x=voltages))

basis.plot(result, shading="gouraud", colorbar=True).show()

current = basis_current.project(
    np.linalg.norm(-basis0.interpolate(sigma) * grad(basis.interpolate(result)), axis=0)
)
basis_current.plot(current, shading="gouraud", colorbar=True, cmap="inferno").show()
