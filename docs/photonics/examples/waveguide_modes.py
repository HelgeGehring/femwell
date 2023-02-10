# ---
# jupyter:
#   jupytext:
#     formats: py:light,md:myst
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# # Modes of a rectangular waveguide

# + tags=["remove-stderr", "hide-input"]
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.affinity
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.mesh import mesh_from_OrderedDict
from femwell.mode_solver import compute_modes, plot_mode

# -

# We describe the geometry using shapely.
# In this case it's simple: we use a shapely.box for the waveguide.
# For the surrounding we buffer the core and clip it to the part below the waveguide for the box.
# The remaining buffer is used as the clad.
# For the core we set the resolution to 30nm and let it fall of over 500nm

# +
wg_width = 2.5
wg_thickness = 0.3
core = shapely.geometry.box(-wg_width / 2, 0, +wg_width / 2, wg_thickness)
env = shapely.affinity.scale(core.buffer(5, resolution=8), xfact=0.5)

polygons = OrderedDict(
    core=core,
    box=clip_by_rect(env, -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(env, -np.inf, 0, np.inf, np.inf),
)

resolutions = dict(core={"resolution": 0.03, "distance": 0.5})

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=10))
mesh.draw().show()
# -

# On this mesh, we define the epsilon. We do this by setting domainwise the epsilon to the squared refractive index.

# +
basis0 = Basis(mesh, ElementTriP0())
epsilon = basis0.zeros()
for subdomain, n in {"core": 1.9963, "box": 1.444, "clad": 1}.items():
    epsilon[basis0.get_dofs(elements=subdomain)] = n**2
basis0.plot(epsilon, colorbar=True).show()
# -

# And now we call `compute_modes` to calculate the modes of the waveguide we set up.
# As modes can have complex fields as soon as the epsilon gets complex, so we get a complex field for each mode.
# Here we show only the real part of the mode.

# +
lams, basis, xs = compute_modes(
    basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1, order=2, radius=15
)
plot_mode(basis, xs[0].real, colorbar=True, direction="x")
plt.show()
