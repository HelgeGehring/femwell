import tempfile
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry
import shapely.affinity
from shapely.ops import clip_by_rect
from skfem import Mesh, Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.mode_solver import compute_modes, plot_mode
from femwell.mesh import mesh_from_OrderedDict

core = shapely.geometry.box(-.5, -.17, .5, .17)
env = shapely.affinity.scale(core.buffer(5, resolution=8), xfact=.5)

polygons = OrderedDict(
    core=core,
    box=clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(env, -np.inf, 0, np.inf, np.inf)
)

resolutions = dict(
    core={"resolution": .03, "distance": .1}
)

with tempfile.TemporaryDirectory() as tmpdirname:
    mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, filename=f'{tmpdirname}/mesh.msh', default_resolution_max=10))

basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros(dtype=complex) + 1
epsilon[basis0.get_dofs(elements='core')] = 1.9963 ** 2
epsilon[basis0.get_dofs(elements='box')] = 1.444 ** 2

lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1, order=2)
plot_mode(basis, xs[0].real, colorbar=True, direction='x')
plt.show()
