from typing import Dict

import numpy as np
from scipy.sparse.linalg import splu
from skfem import Basis, BilinearForm, ElementTriP0, LinearForm, Mesh, asm, enforce
from skfem.helpers import dot

from femwell.thermal import solve_thermal

"""
implemented like in https://www-users.cse.umn.edu/~arnold/8445.f11/notes.pdf page 81 in the middle of the page
"""


def solve_thermal_transient(
    basis0,
    thermal_conductivity_p0,
    thermal_diffusivity_p0,
    specific_conductivity: Dict[str, float],
    current_densities,
    fixed_boundaries,
    dt,
    steps,
):
    basis, temperature = solve_thermal(
        basis0,
        thermal_conductivity_p0,
        specific_conductivity,
        {domain: current(0) for domain, current in current_densities.items()},
        fixed_boundaries=fixed_boundaries,
    )

    @BilinearForm
    def diffusivity_laplace(u, v, w):
        return w["thermal_conductivity"] * dot(u.grad, v.grad)

    @BilinearForm
    def mass(u, v, w):
        return w["thermal_conductivity"] / w["thermal_diffusivity"] * u * v

    L = diffusivity_laplace.assemble(
        basis,
        thermal_conductivity=basis0.interpolate(thermal_conductivity_p0),
    )
    M = mass.assemble(
        basis,
        thermal_diffusivity=basis0.interpolate(thermal_diffusivity_p0),
        thermal_conductivity=basis0.interpolate(thermal_conductivity_p0),
    )

    theta = 0.5  # Crank–Nicolson

    x = basis.zeros()
    for key, value in fixed_boundaries.items():
        x[basis.get_dofs(key)] = value
    L0, M0 = enforce(L, M, D=basis.get_dofs(set(fixed_boundaries.keys())), x=x)
    A = M0 + theta * L0 * dt
    B = M0 - (1 - theta) * L0 * dt

    backsolve = splu(A.T).solve  # .T as splu prefers CSC

    t = 0
    temperatures = [temperature]
    for i in range(steps):
        joule_heating_rhs = basis.zeros()
        for (
            domain,
            current_density,
        ) in current_densities.items():  # sum up the sources for the heating
            current_density_p0 = basis0.zeros()
            current_density_p0[basis0.get_dofs(elements=domain)] = current_density(t)

            @LinearForm
            def joule_heating(v, w):
                return w["current_density"] ** 2 / specific_conductivity[domain] * v

            joule_heating_rhs += asm(
                joule_heating,
                basis,
                current_density=basis0.interpolate(current_density_p0),
            )

        t, temperature = t + dt, backsolve(B @ temperature + joule_heating_rhs * dt)
        temperatures.append(temperature)

    return basis, temperatures


if __name__ == "__main__":
    from collections import OrderedDict

    import matplotlib.pyplot as plt
    from shapely.geometry import LineString, Polygon
    from skfem.io import from_meshio

    from femwell.mesh import mesh_from_OrderedDict

    # Simulating the TiN TOPS heater in https://doi.org/10.1364/OE.27.010456

    w_sim = 8 * 2
    h_clad = 2.8
    h_box = 1
    w_core = 0.5
    h_core = 0.22
    offset_heater = 2.2
    h_heater = 0.14
    w_heater = 2
    h_silicon = 3

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
        heater=Polygon(
            [
                (-w_heater / 2, offset_heater),
                (-w_heater / 2, offset_heater + h_heater),
                (w_heater / 2, offset_heater + h_heater),
                (w_heater / 2, offset_heater),
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
        heater={"resolution": 0.05, "distance": 1},
    )

    mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.3))

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    thermal_conductivity_p0 = basis0.zeros()
    for domain, value in {
        "core": 148,
        "box": 1.38,
        "clad": 1.38,
        "heater": 28,
    }.items():  # , 'silicon': 28
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
    thermal_diffusivity_p0 *= 1e12  # 1e-12 -> conversion from m^2 -> um^2

    dt = 0.1e-5
    steps = 100
    current = (
        lambda t: 0.007 / polygons["heater"].area * ((t < dt * steps / 10) + (t > dt * steps / 2))
    )
    basis, temperatures = solve_thermal_transient(
        basis0,
        thermal_conductivity_p0,
        thermal_diffusivity_p0,
        specific_conductivity={"heater": 2.3e6},
        current_densities={"heater": current},
        fixed_boundaries={"bottom": 0},
        dt=dt,
        steps=steps,
    )

    @LinearForm
    def unit_load(v, w):
        return v

    M = unit_load.assemble(basis)

    times = np.array([dt * i for i in range(steps + 1)])
    plt.xlabel("Time [us]")
    plt.ylabel("Average temperature")
    plt.plot(times * 1e6, M @ np.array(temperatures).T / np.sum(M))
    plt.show()

    # for i in range(0, steps, 10):
    #     fig, ax = plt.subplots(subplot_kw=dict(aspect=1))
    #     for subdomain in mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
    #         mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
    #     basis.plot(temperatures[i], ax=ax, vmin=0, vmax=np.max(temperatures), shading='gouraud').show()

    # Calculate modes

    from tqdm.auto import tqdm

    neffs = []
    for i, temperature in enumerate(tqdm(temperatures)):
        # basis.plot(temperature, vmin=0, vmax=np.max(temperatures))
        # plt.show()

        from femwell.mode_solver import compute_modes, plot_mode

        temperature0 = basis0.project(basis.interpolate(temperature))
        epsilon = basis0.zeros() + (1.444 + 1.00e-5 * temperature0) ** 2
        epsilon[basis0.get_dofs(elements="core")] = (
            3.4777 + 1.86e-4 * temperature0[basis0.get_dofs(elements="core")]
        ) ** 2
        # basis0.plot(epsilon, colorbar=True).show()

        lams, basis_modes, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1)

        # plot_mode(basis_modes, xs[0])
        # plt.show()

        neffs.append(np.real(lams[0]))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Time [us]")
    ax.set_ylabel("Current [mA]")
    ax.plot(times * 1e6, current(times), "b-o")
    ax2 = ax.twinx()
    ax2.set_ylabel("Phase shift")
    ax2.plot(times * 1e6, 2 * np.pi / 1.55 * (neffs - neffs[0]) * 320, "r-o")
    plt.show()
