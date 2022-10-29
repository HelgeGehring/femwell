import tempfile
from tqdm.auto import tqdm

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from skfem import Mesh, Basis, ElementTriP0, ElementVector

from waveguidemodes.mode_solver import compute_modes
from waveguidemodes.waveguide import mesh_waveguide

if __name__ == '__main__':
    # basis0.plot(epsilon, colorbar=True).show()

    num_modes = 10
    widths = np.linspace(.3, 1.5, 100)
    all_lams = np.zeros((widths.shape[0], num_modes))
    all_te_fracs = np.zeros((widths.shape[0], num_modes))
    for i, width in enumerate(tqdm(widths)):
        with tempfile.TemporaryDirectory() as tmpdirname:
            mesh = mesh_waveguide(wsim=2, hclad=.7, hbox=.5, wcore=width, hcore=0.22, filename=tmpdirname + '/mesh.msh')
            mesh = Mesh.load(tmpdirname + '/mesh.msh')

        basis0 = Basis(mesh, ElementTriP0(), intorder=4)
        epsilon = basis0.zeros()
        epsilon[basis0.get_dofs(elements='core')] = 3.4777 ** 2
        epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
        epsilon[basis0.get_dofs(elements='box')] = 1.444 ** 2

        lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=num_modes)
        all_lams[i] = np.real(lams[:num_modes])

        xs = np.real(xs)

        for idx in range(num_modes):
            (et, et_basis), (ez, ez_basis), *_ = basis.split(xs[:, idx])
            plot_basis = et_basis.with_element(ElementVector(ElementTriP0()))
            et_xy = plot_basis.project(et_basis.interpolate(et))
            (et_x, et_x_basis), (et_y, et_y_basis) = plot_basis.split(et_xy)
            all_te_fracs[i, idx] = np.sum(et_x ** 2, axis=0) / np.sum(et_x ** 2 + et_y ** 2, axis=0)

    for lams, te_fracs in zip(all_lams.T, all_te_fracs.T):
        print(te_fracs, cm.seismic(te_fracs))
        plt.plot(widths, lams)
        plt.scatter(widths, lams, c=cm.seismic(2 * (te_fracs - .5)))
    plt.show()
