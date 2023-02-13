"""Waveguide analysis based on https://doi.org/10.1080/02726340290084012."""
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
import scipy.sparse.linalg
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


def compute_modes(
    basis_epsilon_r,
    epsilon_r,
    wavelength,
    mu_r,
    num_modes,
    order=1,
    metallic_boundaries=False,
    radius=np.inf,
    solver="scipy",
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
        epsilon = w.epsilon * (1 + w.x[0] / radius)

        return (
            1 / mu_r * curl(e_t) * curl(v_t)
            - k0**2 * epsilon * dot(e_t, v_t)
            + 1 / mu_r * dot(grad(e_z), v_t)
            + epsilon * inner(e_t, grad(v_z))
            - epsilon * e_z * v_z
        )

    @BilinearForm(dtype=epsilon_r.dtype)
    def bform(e_t, e_z, v_t, v_z, w):
        return -1 / mu_r * dot(e_t, v_t)

    A = aform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))
    B = bform.assemble(basis, epsilon=basis_epsilon_r.interpolate(epsilon_r))

    if metallic_boundaries:
        lams, xs = solve(
            *condense(-A, -B, D=basis.get_dofs()),
            solver=solver(k=num_modes, sigma=k0**2 * np.max(epsilon_r) ** 2),
        )
    else:
        lams, xs = solve(
            -A,
            -B,
            solver=solver(k=num_modes, sigma=k0**2 * np.max(epsilon_r) ** 2),
        )

    idx = np.abs(np.real(lams)).argsort()[::-1]
    lams = lams[idx]
    xs = xs[:, idx]

    xs = xs.T
    xs[:, basis.split_indices()[1]] /= 1j * np.sqrt(
        lams[:, np.newaxis]
    )  # undo the scaling E_3,new = beta * E_3

    for i, lam in enumerate(lams):
        H = calculate_hfield(basis, xs[i], np.sqrt(lam), omega=k0 * scipy.constants.speed_of_light)
        xs[i] /= np.sqrt(calculate_overlap(basis, xs[i], H, basis, xs[i], H))

    return np.sqrt(lams)[:num_modes] / k0, basis, xs[:num_modes]


def calculate_hfield(basis, xs, beta, omega=1):
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


def calculate_overlap(basis_i, E_i, H_i, basis_j, E_j, H_j):
    @Functional
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


def calculate_coupling_coefficient(basis_epsilon, delta_epsilon, basis, E_i, E_j):
    @Functional
    def overlap(w):
        return w["delta_epsilon"] * (
            dot(np.conj(w["E_i"][0]), w["E_j"][0]) + np.conj(w["E_i"][1]) * w["E_j"][1]
        )

    return overlap.assemble(
        basis,
        E_i=basis.interpolate(E_i),
        E_j=basis.interpolate(E_j),
        delta_epsilon=basis_epsilon.interpolate(delta_epsilon),
    )


def confinement_factor(basis_epsilon, epsilon, basis, E):
    @Functional
    def factor(w):
        return (
            np.sqrt(w["epsilon"]) * dot(np.conj(w["E"][0]), w["E"][0])
            + np.conj(w["E"][1]) * w["E"][1]
        )

    return factor.assemble(
        basis,
        E=basis.interpolate(E),
        epsilon=basis_epsilon.interpolate(epsilon),
    )


def calculate_te_frac(basis, x):
    @Functional
    def ex(w):
        return np.abs(w.E[0][0]) ** 2

    @Functional
    def ey(w):
        return np.abs(w.E[0][1]) ** 2

    ex_sum = ex.assemble(basis, E=basis.interpolate(x))
    ey_sum = ey.assemble(basis, E=basis.interpolate(x))

    return ex_sum / (ex_sum + ey_sum)


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
        basis.mesh.draw(ax=ax, boundaries=True, boundaries_only=True)
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

    mesh = Mesh.load("mesh.msh")
    basis = Basis(mesh, ElementTriN2() * ElementTriP2())
    basis0 = basis.with_element(ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    epsilon[basis0.get_dofs(elements="core")] = 3.4777**2
    epsilon[basis0.get_dofs(elements="clad")] = 1.444**2
    epsilon[basis0.get_dofs(elements="box")] = 1.444**2
    # basis0.plot(epsilon, colorbar=True).show()

    lams, basis, xs = compute_modes(
        basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=6, order=2, radius=3
    )

    print(lams)

    plot_mode(basis, np.real(xs[0]))
    plt.show()
    plot_mode(basis, np.imag(xs[0]))
    plt.show()

    xbs = calculate_hfield(basis, xs[0], lams[0] * (2 * np.pi / 1.55))

    plot_mode(basis, np.real(xbs))
    plt.show()
    plot_mode(basis, np.imag(xbs))
    plt.show()

    integrals = np.zeros((len(lams),) * 2, dtype=complex)
    for i in range(len(lams)):
        for j in range(len(lams)):
            E_i = xs[i]
            E_j = xs[j]
            H_i = calculate_hfield(basis, E_i, lams[i] * (2 * np.pi / 1.55))
            H_j = calculate_hfield(basis, E_j, lams[j] * (2 * np.pi / 1.55))
            integrals[i, j] = calculate_overlap(basis, E_i, H_i, basis, E_j, H_j)

    plt.imshow(np.real(integrals))
    plt.colorbar()
    plt.show()
