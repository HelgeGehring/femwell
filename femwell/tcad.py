import matplotlib.pyplot as plt
import numpy as np
import skfem
from scipy.constants import Boltzmann, elementary_charge
from skfem import Basis, BilinearForm, ElementTriP1, LinearForm, condense, solve
from skfem.helpers import grad, inner

# like desicribed in https://www.politesi.polimi.it/bitstream/10589/94468/3/2014_7_Bortolossi.pdf


def solve_coulomb(basis, epsilon_r, fixed_boundaries, phi_k, phi_n, phi_p, doping):
    @BilinearForm
    def coulomb(u, v, w):
        return w.epsilon * inner(grad(u), grad(v))

    def sigma(w):
        return (
            elementary_charge
            * intrinsic_charge
            * (
                np.exp((w.phi_p - w.phi_k) / v_threshold)
                - np.exp((w.phi_k - w.phi_n) / v_threshold)
            )
        )

    @BilinearForm
    def free_charge_contribution(u, v, w):
        return 1 / v_threshold * sigma(w) * inner(u, v)

    @LinearForm
    def fixed_charge_contribution(v, w):
        return (
            w.epsilon * inner(grad(w.phi_k), grad(v))
            + sigma(w) * v
            + elementary_charge * intrinsic_charge * w.doping * v
        )

    A = coulomb.assemble(basis, epsilon=basis.interpolate(epsilon_r))
    B = free_charge_contribution.assemble(
        basis,
        epsilon=basis_epsilon_r.interpolate(epsilon_r),
        phi_k=basis.interpolate(phi_k),
        phi_n=basis.interpolate(phi_n),
        phi_p=basis.interpolate(phi_p),
    )
    C = fixed_charge_contribution.assemble(
        basis,
        epsilon=basis_epsilon_r.interpolate(epsilon_r),
        phi_k=basis.interpolate(phi_k),
        phi_n=basis.interpolate(phi_n),
        phi_p=basis.interpolate(phi_p),
        doping=basis.interpolate(doping),
    )

    u = basis.zeros()
    for key, value in fixed_boundaries.items():
        u[basis.get_dofs(key)] = value

    return basis, solve(
        *condense(A + B, C, x=u, D={key: basis.get_dofs(key) for key in fixed_boundaries})
    )


def solve_continuity_equations(basis, phi_i, p_i_1, n_i_1, pn):
    sign = -1 if pn == "p" else 1
    print(sign)

    @BilinearForm
    def drift_diffusion(u, v, w):
        diffusion_constant = Boltzmann * temperature / elementary_charge * w.mu

        return (
            elementary_charge
            * diffusion_constant
            * np.exp(sign * w.phi_i / v_threshold)
            * inner(grad(u), grad(v))
        )

    def recombination_function(p_i_1, n_i_1):
        return 1 / (
            carrier_lifetime_n * (p_i_1 + intrinsic_charge)
            + carrier_lifetime_p * (n_i_1 + intrinsic_charge)
        )

    @BilinearForm
    def reaction_recombination(u, v, w):
        if pn == "p":
            return (
                elementary_charge
                * w.p_i_1
                / recombination_function(w.p_i_1, w.n_i_1)
                * np.exp(w.phi_i / v_threshold)
                * inner(u, v)
            )
        return (
            elementary_charge
            * w.n_i_1
            / recombination_function(w.p_i_1, w.n_i_1)
            * np.exp(-w.phi_i / v_threshold)
            * inner(u, v)
        )

    @LinearForm
    def force_recombination(v, w):
        return (
            elementary_charge * intrinsic_charge**2 / recombination_function(w.p_i_1, w.n_i_1) * v
        )

    A = drift_diffusion.assemble(basis, phi_i=basis.interpolate(phi_i), mu=400)
    B = reaction_recombination.assemble(
        basis,
        phi_i=basis.interpolate(phi_i),
        p_i_1=basis.interpolate(p_i_1),
        n_i_1=basis.interpolate(n_i_1),
    )
    C = force_recombination.assemble(
        basis, p_i_1=basis.interpolate(p_i_1), n_i_1=basis.interpolate(n_i_1)
    )

    return basis, solve(
        *condense(
            A + B,
            C,
            D=basis.get_dofs(facets="left") + basis.get_dofs(facets="right"),
            x=basis.zeros()
            + basis.project(
                basis.interpolate(doping) * np.exp(-basis.interpolate(phi_i) / v_threshold)
            ),
        )
    )


if __name__ == "__main__":
    k_b = Boltzmann * 1e4  # m^2 -> cm^2
    temperature = 300
    v_threshold = k_b * temperature / elementary_charge
    intrinsic_charge = 10e10
    carrier_lifetime_n = 10e-5
    carrier_lifetime_p = 10e-5

    mesh = skfem.MeshTri().refined(5).scaled(0.05e-4)
    # fig, ax = plt.subplots()
    # mesh.draw(ax=ax).show()

    basis = Basis(mesh, ElementTriP1())
    basis_epsilon_r = basis  # .with_element(ElementTriP0())
    epsilon = basis_epsilon_r.zeros() + 10
    epsilon *= 8.854e-14  # epsilon_0
    # basis.plot(epsilon).show()

    doping = basis.project(lambda x: 2 * (x[0] > 0.025e-4) - 1) * 1e11
    # doping = basis.project(lambda x: x[0]) * 1e17
    fig, ax = plt.subplots()
    basis.plot(doping, ax=ax, shading="gouraud", colorbar=True).show()

    phi_0 = basis.zeros()
    phi_n = basis.zeros()  # basis.project(lambda x: x[0])
    phi_p = basis.zeros()  # basis.project(lambda x: x[0])

    for i in range(10):
        for i in range(10):
            basis_u, d_phi_0 = solve_coulomb(
                basis_epsilon_r,
                epsilon,
                {"left": 1, "right": 0},
                phi_0,
                phi_n,
                phi_p,
                doping,
            )
            phi_0 += d_phi_0 * 1 / 2**i

            fig, ax = plt.subplots()
            ax = basis_u.mesh.draw(ax=ax)
            basis_u.plot(phi_0, ax=ax, shading="gouraud", colorbar=True)
            # basis_vec.plot(-u_grad, ax=ax)
            plt.show()

        basis_n, phi_n_new = solve_continuity_equations(basis_epsilon_r, phi_0, phi_n, phi_p, "n")
        basis_n.plot(phi_n_new, shading="gouraud", colorbar=True).show()

        basis_p, phi_p_new = solve_continuity_equations(basis_epsilon_r, phi_0, phi_n, phi_p, "p")
        basis_p.plot(phi_p_new, shading="gouraud", colorbar=True).show()

        phi_n, phi_p = phi_n_new, phi_p_new
