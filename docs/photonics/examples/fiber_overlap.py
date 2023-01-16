from collections import OrderedDict

from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import shapely.geometry

from skfem import Mesh, Basis, ElementTriP0, ElementTriP1, Functional
from skfem.io import from_meshio

from femwell.mode_solver import compute_modes, plot_mode
from femwell.mesh import mesh_from_OrderedDict
from femwell.fiber import e_field_gaussian


core = shapely.geometry.box(-.1, -.15, .1, .15)

polygons = OrderedDict(
    core=core,
    clad=core.buffer(15, resolution=4)
)

resolutions = dict(
    core={"resolution": .01, "distance": .1}
)

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=10))

basis0 = Basis(mesh, ElementTriP0(), intorder=4)
epsilon = basis0.zeros().astype(complex)
epsilon[basis0.get_dofs(elements='core')] = 1.9963 ** 2
epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
basis0.plot(np.real(epsilon), colorbar=True).show()

lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1)

fig, axs = plot_mode(basis, np.real(xs[0]), direction='x')
plt.show()

mfds = np.linspace(2, 20, 100)
efficiencies = []

for mfd in tqdm(mfds):
    basis_fiber = basis0.with_element(ElementTriP1())
    x_fiber = basis_fiber.project(lambda x: e_field_gaussian(np.sqrt(x[0] ** 2 + x[1] ** 2), 0, mfd / 2, 1, 1.55),
                                    dtype=np.cfloat)


    # basis_fiber.plot(np.real(x_fiber)).show()

    @Functional(dtype=np.complex64)
    def overlap_integral(w):
        return w['E_i'] * np.conj(w['E_j'])

    efficiency = np.abs(overlap_integral.assemble(basis_fiber,
                                                    E_i=basis.interpolate(xs[0])[0][1],
                                                    E_j=basis_fiber.interpolate(x_fiber))
                        /
                        np.sqrt(
                            overlap_integral.assemble(basis_fiber,
                                                        E_i=basis.interpolate(xs[0])[0][1],
                                                        E_j=basis.interpolate(xs[0])[0][1])
                            )
                        )

    efficiencies.append(efficiency)

plt.plot(mfds, efficiencies)
plt.show()
