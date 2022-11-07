from typing import Dict, Iterator, Tuple
from collections import OrderedDict

import numpy as np
from scipy.sparse.linalg import splu
import matplotlib.pyplot as plt
from shapely.geometry import Polygon

from skfem import asm, ElementTriP0, ElementTriP1, BilinearForm, LinearForm, Basis, Mesh, penalize, enforce
from skfem.helpers import dot

from waveguidemodes.mesh import mesh_from_polygons
from waveguidemodes.thermal import solve_thermal

"""
implemented like in https://www-users.cse.umn.edu/~arnold/8445.f11/notes.pdf page 81 in the middle of the page
"""


def solve_thermal_transient(
        basis0,
        thermal_conductivity_p0,
        thermal_diffusivity_p0,
        specific_conductivity: Dict[str, float],
        current_densities,
        dt,
        steps
):
    basis, temperature = solve_thermal(basis0, thermal_conductivity_p0, specific_conductivity,
                                       {domain: current(0) for domain, current in current_densities.items()})

    @BilinearForm
    def diffusivity_laplace(u, v, w):
        return w["thermal_conductivity"] * dot(u.grad, v.grad)

    @BilinearForm
    def mass(u, v, w):
        return w["thermal_conductivity"] / w["thermal_diffusivity"] * u * v

    L = diffusivity_laplace.assemble(
        basis,
        thermal_diffusivity=basis0.interpolate(thermal_diffusivity_p0),
        thermal_conductivity=basis0.interpolate(thermal_conductivity_p0),
    )
    M = mass.assemble(basis,
                      thermal_diffusivity=basis0.interpolate(thermal_diffusivity_p0),
                      thermal_conductivity=basis0.interpolate(thermal_conductivity_p0)
                      )

    theta = 0.5  # Crank–Nicolson

    L0, M0 = penalize(L, M, D=basis.get_dofs(mesh.boundaries["box_None_14"]))
    print(L.shape, L0.shape)
    # L0, M0 = penalize(L, M, D=basis.get_dofs(lambda x: x[1] == np.min(basis.mesh.p[1])))
    A = M0 + theta * L0 * dt
    B = M0 - (1 - theta) * L0 * dt

    backsolve = splu(A.T).solve  # .T as splu prefers CSC

    t = 0
    temperatures = [temperature]
    for i in range(steps):
        joule_heating_rhs = basis.zeros()
        for domain, current_density in current_densities.items():  # sum up the sources for the heating
            current_density_p0 = basis0.zeros()
            current_density_p0[basis0.get_dofs(elements=domain)] = current_density(t)

            @LinearForm
            def joule_heating(v, w):
                return w['current_density'] ** 2 / specific_conductivity[domain]

            joule_heating_rhs += asm(joule_heating, basis,
                                     current_density=basis0.interpolate(current_density_p0)
                                     )
            # basis.plot(joule_heating_rhs).show()

        t, temperature = t + dt, backsolve(B @ temperature + joule_heating_rhs * dt)
        # temperature[basis.get_dofs(mesh.boundaries["box_None_14"])] = 0
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
    h_silicon = 3

    polygons = OrderedDict(
        core=Polygon([
            (-w_core / 2, 0),
            (-w_core / 2, h_core),
            (w_core / 2, h_core),
            (w_core / 2, 0),
        ]),
        heater=Polygon([
            (-w_heater / 2, offset_heater),
            (-w_heater / 2, offset_heater + h_heater),
            (w_heater / 2, offset_heater + h_heater),
            (w_heater / 2, offset_heater),
        ]),
        clad=Polygon([
            (-w_sim / 2, 0),
            (-w_sim / 2, h_clad),
            (w_sim / 2, h_clad),
            (w_sim / 2, 0),
        ]),
        box=Polygon([
            (-w_sim / 2, 0),
            (-w_sim / 2, - h_box),
            (w_sim / 2, - h_box),
            (w_sim / 2, 0),
        ]),
        # silicon=Polygon([
        #    (-w_sim / 2, - h_box - h_silicon),
        #    (-w_sim / 2, - h_box),
        #    (w_sim / 2, - h_box),
        #    (w_sim / 2, - h_box - h_silicon),
        # ]),
    )

    resolutions = dict(
        core={"resolution": 0.05, "distance": 1},
        clad={"resolution": 1, "distance": 1},
        box={"resolution": 1, "distance": 1},
        silicon={"resolution": 1, "distance": 1},
        heater={"resolution": 0.05, "distance": 1}
    )

    mesh_from_polygons(polygons, resolutions, filename='mesh.msh', default_resolution_max=.3)

    mesh = Mesh.load('mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    thermal_conductivity_p0 = basis0.zeros()
    for domain, value in {"core": 148, "box": 1.38, "clad": 1.38, "heater": 28}.items():  # , 'silicon': 28
        thermal_conductivity_p0[basis0.get_dofs(elements=domain)] = value
    thermal_conductivity_p0 *= 1e-12  # 1e-12 -> conversion from 1/m^2 -> 1/um^2

    thermal_diffusivity_p0 = basis0.zeros()
    for domain, value in {
        "heater": 28 / 598 / 5240,
        "box": 1.38 / 709 / 2203,
        "clad": 1.38 / 709 / 2203,
        "core": 148 / 711 / 2330,
        # "silicon": 148 / 711 / 2330,
    }.items():
        thermal_diffusivity_p0[basis0.get_dofs(elements=domain)] = value
    # thermal_diffusivity_p0[:] = 1.38 / 709 / 2203

    thermal_diffusivity_p0 *= 1e12  # 1e-12 -> conversion from m^2 -> um^2

    dt = .1e-6
    steps = 2000
    current = lambda t: 0.007 / polygons['heater'].area * ((t < dt * steps / 10) + (t > dt * steps / 2))
    basis, temperatures = solve_thermal_transient(basis0, thermal_conductivity_p0, thermal_diffusivity_p0,
                                                  specific_conductivity={"heater": 2.3e6},
                                                  current_densities={"heater": current},
                                                  dt=dt,
                                                  steps=steps
                                                  )

    basis.plot(temperatures[0] / temperatures[-1], colorbar=True)
    plt.show()


    @LinearForm
    def unit_load(v, w):
        return v


    M = asm(unit_load, basis)

    print(np.max(temperatures), np.max(temperatures[0]), np.max(temperatures[-1]))
    times = np.array([dt * i for i in range(steps+1)])
    plt.xlabel('Time [us]')
    plt.ylabel('Average temperature')
    plt.plot(times * 1e6, M @ np.array(temperatures).T / np.sum(M))
    plt.show()

    for i in range(0, steps, 100):
        basis.plot(temperatures[i], vmin=0, vmax=np.max(temperatures)).show()

    # Calculate modes
    """
    from tqdm.auto import tqdm

    neffs = []
    for i, temperature in enumerate(tqdm(temperatures)):
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

        # plot_mode(basis_modes, xs[0])
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
    """
