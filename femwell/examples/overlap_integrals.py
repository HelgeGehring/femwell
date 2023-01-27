import tempfile

import matplotlib.pyplot as plt
import numpy as np
from skfem import Basis, ElementTriP0, Mesh
from tqdm.auto import tqdm

from femwell.mode_solver import (
    calculate_hfield,
    calculate_overlap,
    compute_modes,
    plot_mode,
)
from femwell.waveguide import mesh_waveguide

widths = [0.5, 0.8]

modes = []

for width in widths:
    with tempfile.TemporaryDirectory() as tmpdirname:
        mesh = mesh_waveguide(
            wsim=2,
            hclad=0.7,
            hbox=0.5,
            wcore=width,
            hcore=0.22,
            filename=f"{tmpdirname}/mesh.msh",
        )
        mesh = Mesh.load(f"{tmpdirname}/mesh.msh")

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros()
    epsilon[basis0.get_dofs(elements="core")] = 3.4777**2
    epsilon[basis0.get_dofs(elements="clad")] = 1.444**2
    epsilon[basis0.get_dofs(elements="box")] = 1.444**2
    # basis0.plot(epsilon, colorbar=True).show()

    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=10)

    fig, axs = plot_mode(basis, np.real(xs[0]), colorbar=False)
    plt.show()

    modes.append([lams, basis, xs])

for mode_a, mode_b in zip(modes[:-1], modes[1:]):
    lams_a, basis_a, xs_a = mode_a
    lams_b, basis_b, xs_b = mode_b

    integrals = np.zeros((len(lams_a), len(lams_b)), dtype=complex)
    for i in range(len(lams_a)):
        for j in range(len(lams_b)):
            E_i = xs_a[i]
            E_j = xs_b[j]
            H_i = calculate_hfield(basis_a, E_i, lams_a[i] * (2 * np.pi / 1.55))
            H_j = calculate_hfield(basis_b, E_j, lams_b[j] * (2 * np.pi / 1.55))
            integrals[i, j] = calculate_overlap(basis_a, E_i, H_i, basis_b, E_j, H_j)

    plt.imshow(np.abs(integrals))
    plt.colorbar()
    plt.show()

    z = np.abs(integrals) ** 2

    fig, axCenter = plt.subplots(figsize=(8, 8))
    fig.subplots_adjust(0.05, 0.1, 0.95, 0.95)

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(axCenter)
    axvert = divider.append_axes("right", size="30%", pad=0.5)
    axhoriz = divider.append_axes("top", size="20%", pad=0.25)

    axCenter.imshow(z, origin="lower", cmap="jet")
    axhoriz.plot(range(np.shape(z)[1]), np.sum(z, 0))
    axvert.plot(np.sum(z, 1), range(np.shape(z)[0]))

    axhoriz.margins(x=0)
    axvert.margins(y=0)

    plt.show()
