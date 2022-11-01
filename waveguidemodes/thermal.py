from typing import Dict
from collections import OrderedDict

import numpy as np
from shapely.geometry import Polygon

from skfem import asm, ElementTriP0, ElementTriP1, BilinearForm, LinearForm, Basis, solve, condense, Mesh
from skfem.helpers import dot

from waveguidemodes.mesh import mesh_from_polygons


def solve_thermal(
        basis0,
        thermal_conductivity,
        specific_conductivity: Dict[str, float],
        currents,
):
    """Thermal simulation.

    Args:
        Basis: Baisis of the thermal_conductivity
        thermal_conductivity: thermal conductivity in W/m‧K.
        specific_conductivity: specific conductivity in S/m.
        currents: current flowing through the layer in A.

    Returns:
        nothing yet, WIP.
    """

    @BilinearForm
    def conduction(u, v, w):
        return dot(w["thermal_conductivity"] * u.grad, v.grad)

    @LinearForm
    def unit_load(v, _):
        return v

    basis = basis0.with_element(ElementTriP1())
    joule_heating_rhs = basis.zeros()
    for domain, current in currents.items():  # sum up the sources for the heating
        core_basis = Basis(mesh, basis.elem, elements=mesh.subdomains[domain])
        asm_core_unit_load = asm(unit_load, core_basis)
        core_area = np.sum(asm_core_unit_load)
        joule_heating = (current / core_area) ** 2 / specific_conductivity[domain]
        joule_heating_rhs += joule_heating * asm_core_unit_load

    thermal_conductivity_lhs = asm(
        conduction,
        basis,
        thermal_conductivity=basis0.interpolate(thermal_conductivity),
    )

    temperature = solve(
        *condense(
            thermal_conductivity_lhs,
            joule_heating_rhs,
            D=basis.get_dofs(mesh.boundaries["box_None_14"]),
        )
    )

    return basis, temperature


if __name__ == '__main__':
    w_sim = 8 * 2
    h_clad = 2.8
    h_box = 1
    w_core = 0.5
    h_core = 0.22
    offset_heater = 2.2
    h_heater = .14
    w_heater = 2

    polygons = OrderedDict(
        core=Polygon([
            (-w_core / 2, -h_core / 2),
            (-w_core / 2, h_core / 2),
            (w_core / 2, h_core / 2),
            (w_core / 2, -h_core / 2),
        ]),
        heater=Polygon([
            (-w_heater / 2, -h_heater / 2 + offset_heater),
            (-w_heater / 2, h_heater / 2 + offset_heater),
            (w_heater / 2, h_heater / 2 + offset_heater),
            (w_heater / 2, -h_heater / 2 + offset_heater),
        ]),
        clad=Polygon([
            (-w_sim / 2, -h_core / 2),
            (-w_sim / 2, -h_core / 2 + h_clad),
            (w_sim / 2, -h_core / 2 + h_clad),
            (w_sim / 2, -h_core / 2),
        ]),
        box=Polygon([
            (-w_sim / 2, -h_core / 2),
            (-w_sim / 2, -h_core / 2 - h_box),
            (w_sim / 2, -h_core / 2 - h_box),
            (w_sim / 2, -h_core / 2),
        ])
    )

    resolutions = dict(
        core={"resolution": 0.1, "distance": 10},
        clad={"resolution": 0.1, "distance": 10},
        box={"resolution": 0.1, "distance": 10},
        heater={"resolution": 0.01, "distance": 10}
    )

    mesh_from_polygons(polygons, resolutions, filename='mesh.msh', default_resolution_max=.1)

    mesh = Mesh.load('mesh.msh')

    basis0 = Basis(mesh, ElementTriP0())
    thermal_conductivity_p0 = basis0.zeros()
    for domain, value in {"core": 28, "box": 1.38, "clad": 1.38, "heater": 148}.items():
        thermal_conductivity_p0[basis0.get_dofs(elements=domain)] = value
    thermal_conductivity_p0 *= 1e-12  # 1e-12 -> conversion from 1/m^2 -> 1/um^2

    basis, temperature = solve_thermal(basis0, thermal_conductivity_p0,
                                       specific_conductivity={"heater": 2.3e6},
                                       currents={"heater": 0.007})
