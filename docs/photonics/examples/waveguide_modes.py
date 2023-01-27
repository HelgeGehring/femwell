from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.affinity
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0, Mesh
from skfem.io.meshio import from_meshio

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver import compute_modes, plot_mode

wg_width = 0.2
wg_thickness = 0.3
core = shapely.geometry.box(-wg_width / 2, -wg_thickness / 2, +wg_width / 2, +wg_thickness / 2)
env = shapely.affinity.scale(core.buffer(5, resolution=8), xfact=0.5)

polygons = OrderedDict(
    core=core,
    box=clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
)

resolutions = dict(core={"resolution": 0.03, "distance": 0.1})

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=10))

basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros(dtype=complex)
for subdomain, n in {"core": 1.9963, "box": 1.444, "clad": 1}.items():
    epsilon[basis0.get_dofs(elements=subdomain)] = n**2

lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1, order=2)
plot_mode(basis, xs[0].real, colorbar=True, direction="x")
plt.show()
