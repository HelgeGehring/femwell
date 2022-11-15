"""
Adapted from F. Laporte at https://github.com/flaport/meow/blob/main/meow/eme/common.py
"""
import numpy as np
from femwell.mode_solver import calculate_overlap, calculate_scalar_product, calculate_hfield


def compute_interface_s_matrix(
    mode_a, # lams, basis, xs
    mode_b,
):
    lams_a, basis_a, xs_a = mode_a
    lams_b, basis_b, xs_b = mode_b

    products_ab = np.zeros((len(lams_a), len(lams_b)), dtype=complex)
    products_ba = np.zeros((len(lams_b), len(lams_a)), dtype=complex)
    for i in range(len(lams_a)):
        for j in range(len(lams_b)):
            E_i = xs_a[i]
            E_j = xs_b[j]
            H_i = calculate_hfield(basis_a, E_i, -lams_a[i] * (2 * np.pi / 1.55))
            H_j = calculate_hfield(basis_b, E_j, -lams_b[j] * (2 * np.pi / 1.55))
            products_ab[i, j] = np.abs(calculate_scalar_product(basis_a, E_i, basis_b, H_j))
            products_ba[j, i] = np.abs(calculate_scalar_product(basis_b, E_j, basis_a, H_i))

    T_ab = 2*np.linalg.inv(products_ab + products_ba.T)
    R_ab = 0.5*(products_ba.T - products_ab)@T_ab

    T_ba = 2*np.linalg.inv(products_ba + products_ab.T)
    R_ba = 0.5*(products_ab.T - products_ba)@T_ba

    return T_ab, R_ab, T_ba, R_ba


def compute_propagation_s_matrix(modes, length, wavelength):
    lams, basis, xs = modes # lams is neff
    betas = lams * 2 * np.pi / wavelength
    return np.diag(np.exp(2j * np.pi * np.abs(betas) * length))


if __name__ == "__main__":

    import tempfile
    from tqdm.auto import tqdm

    import matplotlib.pyplot as plt
    from matplotlib import cm
    import numpy as np

    from skfem import Mesh, Basis, ElementTriP0

    from femwell.mode_solver import compute_modes, plot_mode, calculate_overlap, calculate_hfield
    from femwell.waveguide import mesh_waveguide

    widths = [.5, .8]

    modes = []

    num_modes = 4
    for width in widths:
        with tempfile.TemporaryDirectory() as tmpdirname:
            mesh = mesh_waveguide(wsim=2, hclad=.7, hbox=.5, wcore=width, hcore=0.22, filename=tmpdirname + '/mesh.msh')
            mesh = Mesh.load(tmpdirname + '/mesh.msh')

        basis0 = Basis(mesh, ElementTriP0(), intorder=4)
        epsilon = basis0.zeros()
        epsilon[basis0.get_dofs(elements='core')] = 3.4777 ** 2
        epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
        epsilon[basis0.get_dofs(elements='box')] = 1.444 ** 2
        #basis0.plot(epsilon, colorbar=True).show()

        lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=num_modes)

        #fig, axs = plot_mode(basis, np.real(xs[0]), colorbar=False)
        # plt.show()

        modes.append([lams, basis, xs])
        print(width, len(lams))

    T_ab, R_ab, T_ba, R_ba = compute_interface_s_matrix(modes[0], modes[1])

    # from pprint import pprint
    # for matrix in [T_ab, R_ab, T_ba, R_ba]:
    #     print(matrix)        

    print(np.shape(compute_propagation_s_matrix(modes[0], 1, 1.55)))
    print(compute_propagation_s_matrix(modes[0], 1, 1.55))