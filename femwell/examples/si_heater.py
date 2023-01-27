from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import LineString, Polygon
from skfem import Basis, ElementTriP0, Mesh

from femwell.mesh import mesh_from_OrderedDict
from femwell.thermal import solve_thermal

# Simulating the TiN TOPS heater in https://doi.org/10.1364/OE.27.010456

w_sim = 8 * 4
h_clad = 2.8
h_box = 2
w_core = 0.5
h_core = 0.22
w_buffer = 0.8
h_heater = h_buffer = 0.09
w_heater = 1
offset_heater = 2.2

polygons = OrderedDict(
    bottom=LineString([(-w_sim / 2, -h_box), (w_sim / 2, -h_box)]),
    core=Polygon(
        [
            (-w_core / 2, 0),
            (-w_core / 2, h_core),
            (w_core / 2, h_core),
            (w_core / 2, 0),
        ]
    ),
    slab_l=Polygon(
        [
            (-w_core / 2 - w_buffer, 0),
            (-w_core / 2 - w_buffer, h_heater),
            (-w_core / 2, h_heater),
            (-w_core / 2, 0),
        ]
    ),
    slab_r=Polygon(
        [
            (+w_core / 2 + w_buffer, 0),
            (+w_core / 2 + w_buffer, h_heater),
            (+w_core / 2, h_heater),
            (+w_core / 2, 0),
        ]
    ),
    heater_l=Polygon(
        [
            (-w_core / 2 - w_buffer - w_heater, 0),
            (-w_core / 2 - w_buffer - w_heater, h_heater),
            (-w_core / 2 - w_buffer, h_heater),
            (-w_core / 2 - w_buffer, 0),
        ]
    ),
    heater_r=Polygon(
        [
            (w_core / 2 + w_buffer + w_heater, 0),
            (w_core / 2 + w_buffer + w_heater, h_heater),
            (w_core / 2 + w_buffer, h_heater),
            (w_core / 2 + w_buffer, 0),
        ]
    ),
    clad=Polygon(
        [
            (-w_sim / 2, 0),
            (-w_sim / 2, h_clad),
            (w_sim / 2, h_clad),
            (w_sim / 2, 0),
        ]
    ),
    box=Polygon(
        [
            (-w_sim / 2, 0),
            (-w_sim / 2, -h_box),
            (w_sim / 2, -h_box),
            (w_sim / 2, 0),
        ]
    ),
)

resolutions = dict(
    core={"resolution": 0.01, "distance": 1},
    clad={"resolution": 0.4, "distance": 1},
    box={"resolution": 0.4, "distance": 1},
    heater_l={"resolution": 0.01, "distance": 1},
    heater_r={"resolution": 0.01, "distance": 1},
)

mesh_from_OrderedDict(polygons, resolutions, filename="mesh.msh", default_resolution_max=0.4)

mesh = Mesh.load("mesh.msh")

basis0 = Basis(mesh, ElementTriP0(), intorder=4)
thermal_conductivity_p0 = basis0.zeros()
for domain, value in {
    "core": 90,
    "box": 1.38,
    "clad": 1.38,
    "slab_l": 55,
    "slab_r": 55,
    "heater_l": 55,
    "heater_r": 55,
}.items():
    thermal_conductivity_p0[basis0.get_dofs(elements=domain)] = value
thermal_conductivity_p0 *= 1e-12  # 1e-12 -> conversion from 1/m^2 -> 1/um^2

power = 25.2e-3
current = np.sqrt(
    power * 1e5 * (polygons["heater_l"].area + polygons["heater_r"].area) * 1e-12 / 320e-6
)
print(current)

basis, temperature = solve_thermal(
    basis0,
    thermal_conductivity_p0,
    specific_conductivity={"heater_l": 1e5, "heater_r": 1e5},
    current_densities={
        "heater_l": current / (polygons["heater_l"].area + polygons["heater_r"].area),
        "heater_r": current / (polygons["heater_l"].area + polygons["heater_r"].area),
    },
    fixed_boundaries={"bottom": 303},
)

fig, ax = plt.subplots(subplot_kw=dict(aspect=1))
for subdomain in mesh.subdomains.keys() - {"gmsh:bounding_entities"}:
    mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
basis.plot(temperature, ax=ax)

from mpl_toolkits.axes_grid1 import make_axes_locatable

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(ax.collections[0], cax=cax)
plt.show()
