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
        dt,
        steps
):
    @BilinearForm
    def conduction(u, v, w):
        return dot(w["thermal_conductivity"] * u.grad, v.grad)

    @LinearForm
    def unit_load(v, _):
        return v

    basis, temperature = solve_thermal(basis0, thermal_conductivity, specific_conductivity,
                                       {domain: current(0) for domain, current in currents.items()})

    basis = basis0.with_element(ElementTriP1())

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

    theta = 0.5  # Crank–Nicolson

    L0, M0 = penalize(L, M, D=basis.get_dofs(mesh.boundaries["box_None_14"]))
    A = M0 + theta * L0 * dt
    B = M0 - (1 - theta) * L0 * dt

    backsolve = splu(A.T).solve  # .T as splu prefers CSC

    t = 0
    temperatures = []
    for i in range(steps):
        joule_heating_rhs = calc_joule_conductivity_rhs(basis, specific_conductivity,
                                                        {domain: current(t) for domain, current in
                                                         currents.items()})
        t, temperature = t + dt, backsolve(B @ temperature + joule_heating_rhs * dt)
        temperatures.append(temperature)

    return basis, temperatures


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
        core={"resolution": 0.05, "distance": 1},
        clad={"resolution": 1, "distance": 1},
        box={"resolution": 1, "distance": 1},
        heater={"resolution": 0.05, "distance": 1}
    )

    mesh_from_polygons(polygons, resolutions, filename='mesh.msh', default_resolution_max=.4)

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

    dt = .4e-6
    steps = 200
    current = lambda t: 0.007 * ((t < dt * steps / 10) + (t > dt * steps / 2))
    basis, temperatures = solve_thermal_transient(basis0, thermal_conductivity_p0,
                                                  specific_conductivity={"heater": 2.3e6},
                                                  currents={"heater": current},
                                                  dt=dt,
                                                  steps=steps
                                                  )

    times = np.array([dt * i for i in range(steps)])
    plt.plot(times * 1e6, np.mean(temperatures, axis=-1))
    plt.show()

    from tqdm.auto import tqdm

    neffs = []
    for temperature in tqdm(temperatures):
        # basis.plot(temperature, vmin=0, vmax=np.max(temperatures))
        # plt.show()

        from waveguidemodes.mode_solver import compute_modes, plot_mode

        temperature0 = basis0.project(basis.interpolate(temperature))
        epsilon = basis0.zeros() + (1.444 + 1.00e-5 * temperature0) ** 2
        epsilon[basis0.get_dofs(elements='core')] = \
            (3.4777 + 1.86e-4 * temperature0[basis0.get_dofs(elements='core')]) ** 2
        # basis0.plot(epsilon, colorbar=True).show()

        lams, basis_modes, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1)

        print(lams)

        # plot_mode(basis, xs[0])
        # plt.show()

        neffs.append(np.real(lams[0]))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time [us]')
    ax.set_ylabel('Current [mA]')
    ax.plot(times * 1e6, current(times), 'b-o')
    ax2 = ax.twinx()
    ax2.set_ylabel('Phase shift')
    ax2.plot(times * 1e6, 2 * np.pi / 1.55 * (neffs - neffs[0]) * 320, 'r-o')
    plt.show()
