from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
from skfem import Basis, ElementDG, ElementTriP1
from skfem.io import from_meshio

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver_2d_periodic import plot_periodic, solve_periodic

height = 5.76 / 2 + 5
a = 0.010
b = 0.78
c = 0.2
slab = 0.920 + 5
pml = 3

wavelength = 1
k0 = 2 * np.pi / wavelength

left = shapely.LineString([(0, y) for y in np.linspace(-height, height, 2)])
right = shapely.LineString([(a, y) for y in np.linspace(-height, height, 2)])
top = shapely.LineString([(x, height) for x in np.linspace(0, a, 2)])
bottom = shapely.LineString([(x, -height) for x in np.linspace(0, a, 2)])

box = shapely.box(0, -height, a, height)
structure = shapely.box(0, -b / 2, a, b / 2)
structure1 = shapely.box(0, height - slab, a, height)
structure2 = shapely.box(0, -height + slab, a, -height)

resolutions = {
    "structure": {"resolution": 0.1, "distance": 0.1},
    "hole": {"resolution": 0.1, "distance": 0.1},
}

mesh = from_meshio(
    mesh_from_OrderedDict(
        OrderedDict(
            left=left,
            right=right,
            top=top,
            bottom=bottom,
            structure=structure,
            structure1=structure1,
            structure2=structure2,
            box=box,
        ),
        resolutions=resolutions,
        filename="mesh.msh",
        default_resolution_max=0.05,
    )
)

basis_vec = Basis(mesh, ElementTriP1() * ElementTriP1())
basis_epsilon_r = basis_vec.with_element(ElementDG(ElementTriP1()))

epsilon_r = basis_epsilon_r.zeros(dtype=np.complex64) + 1.45
epsilon_r[basis_epsilon_r.get_dofs(elements="box")] = 1.39
epsilon_r **= 2
basis_epsilon_r.plot(np.real(epsilon_r), aspect=a, colorbar=True).show()

epsilon_r += basis_epsilon_r.project(
    lambda x: (0.5j) * (np.clip(np.abs(x[1]) - height + pml, 0, np.inf) / pml) ** 2,
    dtype=np.complex64,
)
basis_epsilon_r.plot(np.imag(epsilon_r), aspect=a, colorbar=True).show()

ks, basis_phi, phis = solve_periodic(basis_epsilon_r, epsilon_r, k0)

for i, k in enumerate(ks):
    if 0 < np.imag(k) < 0.002:
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))
        plt.title(f"{k}")
        plot_periodic(k, a, basis_phi, phis[..., i], 100, ax)
        plt.show()
