from typing import Dict, Iterator, Tuple
from collections import OrderedDict

import numpy as np
from scipy.sparse.linalg import splu
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from shapely.geometry import Polygon

from skfem import asm, ElementTriP0, ElementTriP1, BilinearForm, LinearForm, Basis, Mesh, penalize
from skfem.helpers import dot
from skfem.visuals.matplotlib import draw, plot

from waveguidemodes.mesh import mesh_from_polygons
from waveguidemodes.thermal import solve_thermal, calc_joule_conductivity_rhs


def solve_thermal_transient(
        basis0,
        thermal_conductivity,
        specific_conductivity: Dict[str, float],
        currents,
):
    @BilinearForm
    def conduction(u, v, w):
        return dot(w["thermal_conductivity"] * u.grad, v.grad)

    @LinearForm
    def unit_load(v, _):
        return v

    basis, temperature = solve_thermal(basis0, thermal_conductivity, specific_conductivity, currents)

    basis = basis0.with_element(ElementTriP1())
    joule_heating_rhs = calc_joule_conductivity_rhs(basis, specific_conductivity, currents)

    @BilinearForm
    def diffusivity_laplace(u, v, w):
        return dot(u.grad * w["thermal_conductivity"], v.grad)

    @BilinearForm
    def mass(u, v, w):
        return u * v / (w["thermal_diffusivity"] / w["thermal_conductivity"])

    L = asm(
        diffusivity_laplace,
        basis,
        thermal_diffusivity=basis0.interpolate(thermal_diffusivity_p0),
        thermal_conductivity=basis0.interpolate(thermal_conductivity_p0),
    )
    M = asm(
        mass,
        basis,
        thermal_diffusivity=basis0.interpolate(thermal_diffusivity_p0),
        thermal_conductivity=basis0.interpolate(thermal_conductivity_p0),
    )

    dt = 0.1e-6
    theta = 0.5  # Crank–Nicolson
    steps = 200

    L0, M0 = penalize(L, M, D=basis.get_dofs(mesh.boundaries["box_None_14"]))
    A = M0 + theta * L0 * dt
    B = M0 - (1 - theta) * L0 * dt

    backsolve = splu(A.T).solve  # .T as splu prefers CSC

    def evolve(
            t: float, u: np.ndarray, heating: np.ndarray
    ) -> Iterator[Tuple[float, np.ndarray]]:
        i = 0
        while True:
            t_temperature[i] = t, np.mean(u)
            i += 1
            yield t, u
            t, u = t + dt, backsolve(B @ u + heating * dt)

    ax = draw(mesh, boundaries_only=True)
    ax.set_axis_on()
    ax = plot(mesh, temperature, ax=ax, shading="gouraud")
    title = ax.set_title("t = 0.00")
    field = ax.get_children()[1]  # vertex-based temperature-colour
    fig = ax.get_figure()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(field, cax=cax)

    def update(event):
        t, u = event
        title.set_text(f"$t$ = {t * 1e6:.2f}us")
        field.set_array(u)

    t_temperature = np.zeros((steps + 2, 2))
    animation = FuncAnimation(
        fig,
        update,
        evolve(0.0, temperature * 0.01, joule_heating_rhs),
        repeat=False,
        interval=30,
        save_count=steps,
    )
    animation.save("heater_up.gif", "imagemagick")
    t_temperature_up = t_temperature

    t_temperature = np.zeros((steps + 2, 2))
    animation = FuncAnimation(
        fig,
        update,
        evolve(0.0, temperature, 0),
        repeat=False,
        interval=30,
        save_count=steps,
    )
    animation.save("heater_down.gif", "imagemagick")
    t_temperature_down = t_temperature

    plt.figure()
    plt.plot(t_temperature[:-1, 0] * 1e6, t_temperature_up[:-1, 1])
    plt.plot(t_temperature[:-1, 0] * 1e6, t_temperature_down[:-1, 1])
    plt.plot(
        t_temperature[:-1, 0] * 1e6, t_temperature[:-1, 1] * 0 + np.mean(temperature)
    )
    plt.xlabel("Time [us]")
    plt.ylabel("Average temperature offset [T]")
    plt.savefig("heater.svg", bbox_inches="tight")
    plt.show()


if __name__ == '__main__':
    # Simulating the TiN TOPS heater in https://doi.org/10.1364/OE.27.010456

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
        core={"resolution": 0.02, "distance": 1},
        clad={"resolution": 0.4, "distance": 1},
        box={"resolution": 0.4, "distance": 1},
        heater={"resolution": 0.05, "distance": 1}
    )

    mesh_from_polygons(polygons, resolutions, filename='mesh.msh', default_resolution_max=.1)

    mesh = Mesh.load('mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    thermal_conductivity_p0 = basis0.zeros()
    for domain, value in {"core": 28, "box": 1.38, "clad": 1.38, "heater": 148}.items():
        thermal_conductivity_p0[basis0.get_dofs(elements=domain)] = value
    thermal_conductivity_p0 *= 1e-12  # 1e-12 -> conversion from 1/m^2 -> 1/um^2

    thermal_diffusivity_p0 = basis0.zeros()
    for domain, value in {
        "heater": 28 / 598 / 5240,
        "box": 1.38 / 709 / 2203,
        "clad": 1.38 / 709 / 2203,
        "core": 148 / 711 / 2330,
    }.items():
        thermal_diffusivity_p0[basis0.get_dofs(elements=domain)] = value

    thermal_diffusivity_p0 *= 1e12  # 1e-12 -> conversion from m^2 -> um^2

    solve_thermal_transient(basis0, thermal_conductivity_p0,
                            specific_conductivity={"heater": 2.3e6},
                            currents={"heater": 0.007})
