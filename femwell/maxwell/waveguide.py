"""Waveguide analysis based on https://doi.org/10.1080/02726340290084012."""
from dataclasses import dataclass
from functools import cached_property
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import scipy.sparse.linalg
from matplotlib.axes import Axes
from matplotlib.figure import Figure
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
    InteriorFacetBasis,
    LinearForm,
    Mesh,
    condense,
    solve,
)
from skfem.helpers import cross, curl, dot, grad, inner
from skfem.utils import solver_eigen_scipy


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
        @Functional(dtype=complex)
        def overlap(w):
            return w["delta_epsilon"] * (
                dot(np.conj(w["E_i"][0]), w["E_j"][0]) + np.conj(w["E_i"][1]) * w["E_j"][1]
            )

        return overlap.assemble(
            self.basis,
            E_i=self.basis.interpolate(self.E),
            E_j=self.basis.interpolate(mode.E),
            delta_epsilon=self.basis_epsilon_r.interpolate(delta_epsilon),
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
        @Functional
        def factor(w):
            return np.sqrt(w["epsilon"]) * (
                dot(np.conj(w["E"][0]), w["E"][0]) + np.conj(w["E"][1]) * w["E"][1]
            )

        basis = self.basis.with_elements(elements)
        basis_epsilon_r = self.basis_epsilon_r.with_elements(elements)
        return (
            speed_of_light
            * epsilon_0
            * factor.assemble(
                basis,
                E=basis.interpolate(self.E),
                epsilon=basis_epsilon_r.interpolate(self.epsilon_r),
            )
        )

    def calculate_pertubated_neff(self, delta_epsilon):
        return (
            self.n_eff
            + self.calculate_coupling_coefficient(self, delta_epsilon)
            * scipy.constants.epsilon_0
            * scipy.constants.speed_of_light
            * 0.5
        )

    def calculate_intensity(self) -> Tuple[NDArray, Basis]:
        """Calculates the intensity of a mode.

        The intensity is calculated from the cross-product between the electric and magnetic field, as
        described in https://doi.org/10.1364/OE.16.016659.

        The calculation is performed as follows:
        1) The electric and magnetic fields are interpolated on the quadrature points with the simulation basis;
        2) The intensity is calculated directly on the quadrature points;
        3) The intensity is projected on a new discontinuous, piecewise linear triangular basis.

        Returns:
            Basis, NDArray: Plot-ready basis and intensity array
        """
        (Ex, Ey), _ = self.basis.interpolate(self.E)
        (Hx, Hy), _ = self.basis.interpolate(np.conj(self.H))
        intensity = 0.5 * np.real(Ex * Hy - Ey * Hx)
        basis2 = self.basis.with_element(ElementDG(ElementTriP1()))
        intensity2 = basis2.project(intensity)

        return basis2, intensity2

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

    def plot_intensity(
        self,
        ax: Axes = None,
        colorbar: bool = True,
        normalize: bool = True,
    ) -> Tuple[Figure, Axes]:
        """Plots the intensity of a mode as outlined in `calculate_intensity`.

        Args:
            ax (Axes, optional): Axes onto which the plot is drawn. Defaults to None.
            colorbar (bool, optional): Adds a colorbar to the plot. Defaults to True.
            normalize (bool, optional): Normalizes the intensity by its maximum value. Defaults to True.

        Returns:
            Tuple[Figure, Axes]: Figure and axes of the plot.
        """
        intensity_basis, intensity = self.calculate_intensity()
        if normalize:
            intensity = intensity / intensity.max()

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        for subdomain in self.basis.mesh.subdomains.keys() - {"gmsh:bounding_entities"}:
            self.basis.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True, color="w")
        intensity_basis.plot(intensity, ax=ax, cmap="inferno")

        if colorbar:
            plt.colorbar(ax.collections[-1])

        return fig, ax


@dataclass(frozen=True)
class Modes:
    modes: List

    def __getitem__(self, idx) -> Mode:
        return self.modes[idx]

    def __len__(self) -> int:
        return len(self.modes)

    def __repr__(self) -> str:
        modes = "\n\t" + "\n\t".join(repr(mode) for mode in self.modes) + "\n"
        return f"{self.__class__.__name__}(modes=({modes}))"

    def sorted(self, key):
        return Modes(modes=sorted(self.modes, key=key))

    @property
    def n_effs(self):
        return np.array([mode.n_eff for mode in self.modes])


def compute_modes(
    basis_epsilon_r,
    epsilon_r,
    wavelength,
    *,
    mu_r=1,
    num_modes=1,
    order=1,
    metallic_boundaries=False,
    radius=np.inf,
    n_guess=None,
    solver="scipy",
) -> Modes:
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
        sigma = sigma = k0**2 * np.max(epsilon_r) * 1.1

    if metallic_boundaries:
        lams, xs = solve(
            *condense(
                -A,
                -B,
                D=basis.get_dofs(None if metallic_boundaries is True else metallic_boundaries),
                x=basis.zeros(dtype=complex),
            ),
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


def calculate_hfield(basis, xs, beta, omega):
    @BilinearForm(dtype=np.complex64)
    def aform(e_t, e_z, v_t, v_z, w):
        return (
            (-1j * beta * e_t[1] + e_z.grad[1]) * v_t[0]
            + (1j * beta * e_t[0] - e_z.grad[0]) * v_t[1]
            + e_t.curl * v_z
        )

    @BilinearForm(dtype=np.complex64)
    def bform(e_t, e_z, v_t, v_z, w):
        return dot(e_t, v_t) + e_z * v_z

    return (
        scipy.sparse.linalg.spsolve(
            bform.assemble(basis), aform.assemble(basis) @ xs.astype(complex)
        )
        * -1j
        / scipy.constants.mu_0
        / omega
    )


def calculate_energy_current_density(basis, xs):
    basis_energy = basis.with_element(ElementTriP0())

    @LinearForm(dtype=complex)
    def aform(v, w):
        e_t, e_z = w["e"]
        return abs(e_t[0]) ** 2 * v + abs(e_t[1]) ** 2 * v + abs(e_z) * v

    a_operator = aform.assemble(basis_energy, e=basis.interpolate(xs))

    @BilinearForm(dtype=complex)
    def bform(e, v, w):
        return e * v

    b_operator = bform.assemble(basis_energy)

    return basis_energy, scipy.sparse.linalg.spsolve(b_operator, a_operator)


def calculate_overlap(
    basis_i: Basis,
    E_i: np.ndarray,
    H_i: np.ndarray,
    basis_j: Basis,
    E_j: np.ndarray,
    H_j: np.ndarray,
) -> np.complex64:
    """Calculates the fully vectorial overlap between two modes.

    If the modes do not share the basis, interpolation is performed automatically.

    Args:
        basis_i (Basis): Basis of the first mode
        E_i (np.ndarray): Electric field of the first mode
        H_i (np.ndarray): Magnetic field of the first mode
        basis_j (Basis): Basis of the second mode
        E_j (np.ndarray): Electric field of the first mode
        H_j (np.ndarray): Magnetic field of the first mode

    Returns:
        np.complex64: Complex overlap between the two modes
    """

    @Functional(dtype=np.complex64)
    def overlap(w):
        return cross(np.conj(w["E_i"][0]), w["H_j"][0]) + cross(w["E_j"][0], np.conj(w["H_i"][0]))

    if basis_i == basis_j:
        return 0.5 * overlap.assemble(
            basis_i,
            E_i=basis_i.interpolate(E_i),
            H_i=basis_i.interpolate(H_i),
            E_j=basis_j.interpolate(E_j),
            H_j=basis_j.interpolate(H_j),
        )
    basis_j_fix = basis_j.with_element(ElementVector(ElementTriP1()))

    (et, et_basis), (ez, ez_basis) = basis_j.split(E_j)
    E_j = basis_j_fix.project(et_basis.interpolate(et), dtype=np.cfloat)
    (et_x, et_x_basis), (et_y, et_y_basis) = basis_j_fix.split(E_j)

    (et, et_basis), (ez, ez_basis) = basis_j.split(H_j)
    H_j = basis_j_fix.project(et_basis.interpolate(et), dtype=np.cfloat)
    (ht_x, ht_x_basis), (ht_y, ht_y_basis) = basis_j_fix.split(H_j)

    @Functional(dtype=np.complex64)
    def overlap(w):
        return cross(
            np.conj(w["E_i"][0]),
            np.array((ht_x_basis.interpolator(ht_x)(w.x), ht_y_basis.interpolator(ht_y)(w.x))),
        ) + cross(
            np.array((et_x_basis.interpolator(et_x)(w.x), et_y_basis.interpolator(et_y)(w.x))),
            np.conj(w["H_i"][0]),
        )

    return 0.5 * overlap.assemble(
        basis_i, E_i=basis_i.interpolate(E_i), H_i=basis_i.interpolate(H_i)
    )


def calculate_scalar_product(basis_i, E_i, basis_j, H_j):
    @Functional
    def overlap(w):
        return cross(np.conj(w["E_i"][0]), w["H_j"][0])

    if basis_i == basis_j:
        return overlap.assemble(basis_i, E_i=basis_i.interpolate(E_i), H_j=basis_j.interpolate(H_j))
    basis_j_fix = basis_j.with_element(ElementVector(ElementTriP1()))

    (et, et_basis), (ez, ez_basis) = basis_j.split(H_j)
    H_j = basis_j_fix.project(et_basis.interpolate(et), dtype=np.cfloat)
    (ht_x, ht_x_basis), (ht_y, ht_y_basis) = basis_j_fix.split(H_j)

    @Functional(dtype=np.complex64)
    def overlap(w):
        return cross(
            np.conj(w["E_i"][0]),
            np.array((ht_x_basis.interpolator(ht_x)(w.x), ht_y_basis.interpolator(ht_y)(w.x))),
        )

    return overlap.assemble(basis_i, E_i=basis_i.interpolate(E_i))


def plot_mode(basis, mode, plot_vectors=False, colorbar=True, title="E", direction="y"):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    (et, et_basis), (ez, ez_basis) = basis.split(mode)

    if plot_vectors:
        rc = (2, 1) if direction != "x" else (1, 2)
        fig, axs = plt.subplots(*rc, subplot_kw=dict(aspect=1))
        for ax in axs:
            basis.mesh.draw(ax=ax, boundaries=True, boundaries_only=True)
            for subdomain in basis.mesh.subdomains.keys() - {"gmsh:bounding_entities"}:
                basis.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
        et_basis.plot(et, ax=axs[0])
        ez_basis.plot(ez, ax=axs[1])

        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(axs[1].collections[0], cax=cax)
        return fig, axs

    plot_basis = et_basis.with_element(ElementVector(ElementDG(ElementTriP1())))
    et_xy = plot_basis.project(et_basis.interpolate(et))
    (et_x, et_x_basis), (et_y, et_y_basis) = plot_basis.split(et_xy)

    rc = (3, 1) if direction != "x" else (1, 3)
    fig, axs = plt.subplots(*rc, subplot_kw=dict(aspect=1))
    for ax in axs:
        basis.mesh.draw(ax=ax, boundaries_only=True)
        for subdomain in basis.mesh.subdomains.keys() - {"gmsh:bounding_entities"}:
            basis.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)

    for ax, component in zip(axs, "xyz"):
        ax.set_title(f"${title}_{component}$")

    maxabs = max(np.max(np.abs(data.value)) for data in basis.interpolate(mode))
    vmin = -maxabs if colorbar == "same" else None
    vmax = maxabs if colorbar == "same" else None

    et_x_basis.plot(et_x, shading="gouraud", ax=axs[0], vmin=vmin, vmax=vmax)
    et_y_basis.plot(et_y, shading="gouraud", ax=axs[1], vmin=vmin, vmax=vmax)
    ez_basis.plot(ez, shading="gouraud", ax=axs[2], vmin=vmin, vmax=vmax)

    if colorbar:
        if colorbar == "same":
            plt.colorbar(axs[0].collections[-1], ax=axs.ravel().tolist())
        else:
            for ax in axs:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(ax.collections[-1], cax=cax)
            plt.tight_layout()

    return fig, axs


def eval_error_estimator(basis, u):
    @Functional
    def interior_residual(w):
        h = w.h
        x, y = w.x
        return h**2  # * load_func(x, y) ** 2

    eta_K = interior_residual.elemental(basis, w=basis.interpolate(u))

    # facet jump
    fbasis = [InteriorFacetBasis(basis.mesh, basis.elem, side=i) for i in [0, 1]]
    w = {"u" + str(i + 1): fbasis[i].interpolate(u) for i in [0, 1]}

    @Functional
    def edge_jump(w):
        return w.h * (
            np.abs(dot(grad(w["u1"][1]) - grad(w["u2"][1]), w.n)) ** 2
            + np.abs(dot(w["u1"][0] - w["u2"][0], w.n)) ** 2
        )

    tmp = np.zeros(basis.mesh.facets.shape[1])
    tmp[fbasis[0].find] = edge_jump.elemental(fbasis[0], **w)
    eta_E = np.sum(0.5 * tmp[basis.mesh.t2f], axis=0)

    return eta_K + eta_E


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
