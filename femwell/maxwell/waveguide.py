"""Waveguide analysis based on https://doi.org/10.1080/02726340290084012."""
from dataclasses import dataclass
from functools import cached_property
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import scipy.sparse.linalg
from numpy.typing import NDArray
from scipy.constants import epsilon_0, speed_of_light
from skfem import (
    Basis,
    BilinearForm,
    ElementDG,
    ElementTriN1,
    ElementTriN2,
    ElementTriP0,
    ElementTriP1,
    ElementTriP2,
    ElementVector,
    Functional,
    LinearForm,
    Mesh,
    condense,
    solve,
)
from skfem.helpers import cross, curl, dot, grad, inner
from skfem.utils import solver_eigen_scipy

from femwell.mode_solver import (
    calculate_coupling_coefficient,
    calculate_energy_current_density,
    calculate_hfield,
    calculate_overlap,
    calculate_scalar_product,
    confinement_factor,
    plot_mode,
)


@dataclass(frozen=True)
class Modes:
    modes: List

    def __getitem__(self, idx):
        return self.modes[idx]

    def __len__(self):
        return len(self.modes)

    def __repr__(self) -> str:
        modes = "\n\t" + "\n\t".join(repr(mode) for mode in self.modes) + "\n"
        return f"{self.__class__.__name__}(modes=({modes}))"

    def sorted(self, key):
        return Modes(modes=sorted(self.modes, key=key))


@dataclass(frozen=True)
class Mode:
    frequency: float
    """Frequency of the light"""
    k: float
    """Propagation constant of the mode"""
    basis_epsilon_r: Basis
    """Basis used for epsilon_r"""
    epsilon_r: NDArray
    """Epsilon_r with which the mode was calculated"""
    basis: Basis
    """Basis on which the mode was calculated and E/H are defined"""
    E: NDArray
    """Electric field of the mode"""
    H: NDArray
    """Magnetic field of the mode"""

    @property
    def omega(self):
        """Angular frequency of the light"""
        return 2 * np.pi * self.frequency

    @property
    def k0(self):
        """Vacuum propagation constant of the light"""
        return self.omega / speed_of_light

    @property
    def wavelength(self):
        """Vacuum wavelength of the light"""
        return speed_of_light / self.frequency

    @property
    def n_eff(self):
        """Effective refractive index of the mode"""
        return self.k / self.k0

    @cached_property
    def te_fraction(self):
        """TE-fraction of the mode"""

        @Functional
        def ex(w):
            return np.abs(w.E[0][0]) ** 2

        @Functional
        def ey(w):
            return np.abs(w.E[0][1]) ** 2

        ex_sum = ex.assemble(self.basis, E=self.basis.interpolate(self.E))
        ey_sum = ey.assemble(self.basis, E=self.basis.interpolate(self.E))

        return ex_sum / (ex_sum + ey_sum)

    @cached_property
    def tm_fraction(self):
        """TM-fraction of the mode"""

        return 1 - self.te_fraction

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(k: {self.k}, n_eff:{self.n_eff})"

    def calculate_overlap(self, mode, elements=None):
        return calculate_overlap(self.basis, self.E, self.H, mode.basis, mode.E, mode.H)

    def calculate_coupling_coefficient(self, mode, delta_epsilon):
        return calculate_coupling_coefficient(
            self.basis_epsilon_r, delta_epsilon, self.basis, self.E, mode.E
        )

    def calculate_propagation_loss(self, distance):
        return -20 / np.log(10) * self.k0 * np.imag(self.n_eff) * distance

    def calculate_power(self, elements=None):
        if not elements:
            basis = self.basis
        else:
            basis = self.basis.with_elements(elements)
        return calculate_overlap(basis, self.E, self.H, basis, self.E, self.H)

    def calculate_confinement_factor(self, elements):
        return confinement_factor(
            self.basis_epsilon_r.with_elements(elements),
            self.epsilon_r,
            self.basis.with_elements(elements),
            self.E,
        )

    def plot(self, field, plot_vectors=False, colorbar=True, direction="y", title="E"):
        return plot_mode(
            self.basis,
            field,
            plot_vectors=plot_vectors,
            colorbar=colorbar,
            title=title,
            direction=direction,
        )

    def show(self, field, **kwargs):
        self.plot(field=field, **kwargs)
        plt.show()


def compute_modes(
    basis_epsilon_r,
    epsilon_r,
    wavelength,
    mu_r=1,
    num_modes=1,
    order=1,
    metallic_boundaries=False,
    radius=np.inf,
    n_guess=None,
    solver="slepc",
):
    if solver == "scipy":
        solver = solver_eigen_scipy
    elif solver == "slepc":
        from femwell.solver import solver_eigen_slepc

        solver = solver_eigen_slepc
    else:
        raise ValueError("`solver` must either be `scipy` or `slepc`")

    k0 = 2 * np.pi / wavelength

    if order == 1:
        element = ElementTriN1() * ElementTriP1()
    elif order == 2:
        element = ElementTriN2() * ElementTriP2()
    else:
        raise AssertionError("Only order 1 and 2 implemented by now.")

    basis = basis_epsilon_r.with_element(element)
    basis_epsilon_r = basis.with_element(basis_epsilon_r.elem)  # adjust quadrature

    @BilinearForm(dtype=epsilon_r.dtype)
    def aform(e_t, e_z, v_t, v_z, w):
        epsilon = w.epsilon * (1 + w.x[0] / radius) ** 2

        return (
            1 / mu_r * curl(e_t) * curl(v_t) / k0**2
            - epsilon * dot(e_t, v_t)
            + 1 / mu_r * dot(grad(e_z), v_t)
            + epsilon * inner(e_t, grad(v_z))
            - epsilon * e_z * v_z * k0**2
        )

    @BilinearForm(dtype=epsilon_r.dtype)
    def bform(e_t, e_z, v_t, v_z, w):
        return -1 / mu_r * dot(e_t, v_t) / k0**2

    A = aform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))
    B = bform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))

    if n_guess:
        sigma = sigma = k0**2 * n_guess**2
    else:
        sigma = sigma = k0**2 * np.max(epsilon_r) ** 2

    if metallic_boundaries:
        lams, xs = solve(
            *condense(-A, -B, D=basis.get_dofs(), x=basis.zeros(dtype=complex)),
            solver=solver(k=num_modes, sigma=sigma),
        )
    else:
        lams, xs = solve(
            -A,
            -B,
            solver=solver(k=num_modes, sigma=sigma),
        )

    xs[basis.split_indices()[1], :] /= 1j * np.sqrt(
        lams[np.newaxis, :] / k0**4
    )  # undo the scaling E_3,new = beta * E_3

    hs = []
    for i, lam in enumerate(lams):
        H = calculate_hfield(
            basis, xs[:, i], np.sqrt(lam), omega=k0 * scipy.constants.speed_of_light
        )
        power = calculate_overlap(basis, xs[:, i], H, basis, xs[:, i], H)
        xs[:, i] /= np.sqrt(power)
        H /= np.sqrt(power)
        hs.append(H)

    return Modes(
        modes=[
            Mode(
                frequency=speed_of_light / wavelength,
                k=np.sqrt(lams[i]),
                basis_epsilon_r=basis_epsilon_r,
                epsilon_r=epsilon_r,
                basis=basis,
                E=xs[:, i],
                H=hs[i],
            )
            for i in range(num_modes)
        ]
    )


if __name__ == "__main__":
    from collections import OrderedDict

    from shapely.geometry import Polygon

    from femwell.mesh import mesh_from_OrderedDict

    x_min = 0
    w_sim = 3
    h_clad = 0.7
    h_box = 0.5
    w_core = 1
    h_core = 0.22
    offset_heater = 2.2
    h_heater = 0.14
    w_heater = 2

    polygons = OrderedDict(
        core=Polygon(
            [
                (x_min - w_core / 2, 0),
                (x_min - w_core / 2, h_core),
                (x_min + w_core / 2, h_core),
                (x_min + w_core / 2, 0),
            ]
        ),
        clad=Polygon(
            [
                (x_min - w_sim / 2, 0),
                (x_min - w_sim / 2, h_clad),
                (x_min + w_sim / 2, h_clad),
                (x_min + w_sim / 2, 0),
            ]
        ),
        box=Polygon(
            [
                (x_min - w_sim / 2, 0),
                (x_min - w_sim / 2, -h_box),
                (x_min + w_sim / 2, -h_box),
                (x_min + w_sim / 2, 0),
            ]
        ),
    )

    resolutions = dict(core={"resolution": 0.05, "distance": 1})

    mesh_from_OrderedDict(polygons, resolutions, filename="mesh.msh", default_resolution_max=0.2)

    scale = 1e0
    mesh = Mesh.load("mesh.msh").scaled(scale)
    basis = Basis(mesh, ElementTriN2() * ElementTriP2())
    basis0 = basis.with_element(ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    epsilon[basis0.get_dofs(elements="core")] = 3.4777**2
    epsilon[basis0.get_dofs(elements="clad")] = 1.444**2
    epsilon[basis0.get_dofs(elements="box")] = 1.444**2
    # basis0.plot(epsilon, colorbar=True).show()

    modes = compute_modes(
        basis0,
        epsilon,
        wavelength=1.55 * scale,
        mu_r=1,
        num_modes=6,
        order=2,
        radius=3 * scale,
    )
    print(modes)
    print(modes[0].te_fraction)

    modes[0].show(np.real(modes[0].E))
    modes[0].show(np.imag(modes[0].E))

    modes[0].show(np.real(modes[0].H))
    modes[0].show(np.imag(modes[0].H))

    integrals = np.zeros((len(modes),) * 2, dtype=complex)

    for i in range(len(modes)):
        for j in range(len(modes)):
            integrals[i, j] = modes[i].calculate_overlap(modes[j])

    plt.imshow(np.real(integrals))
    plt.colorbar()
    plt.show()

    # Create basis to select a certain simulation extent
    def sel_fun(x):
        return (x[0] < 0) * (x[0] > -1) * (x[1] > 0) * (x[1] < 0.5)

    print(modes.sorted(lambda mode: mode.calculate_power(elements=sel_fun)))
