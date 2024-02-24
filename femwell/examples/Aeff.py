from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import box
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from femwell.maxwell.waveguide import compute_modes
from femwell.maxwell.waveguide import calculate_sfwm_Aeff
from femwell.mesh import mesh_from_OrderedDict
from scipy.constants import c, epsilon_0


#Dispersion relations of materials
#Core
def n_X(wavelength):
    x = wavelength
    return (1 
            + 2.19244563 / (1 - (0.20746607 / x) ** 2)
            + 0.13435116 / (1 - (0.3985835 / x) ** 2)
            + 2.20997784 / (1 - (0.20747044 / x) ** 2)
    ) ** 0.5

# Box
def n_silicon_dioxide(wavelength):
    x = wavelength
    return (
        1
        + 0.6961663 / (1 - (0.0684043 / x) ** 2)
        + 0.4079426 / (1 - (0.1162414 / x) ** 2)
        + 0.8974794 / (1 - (9.896161 / x) ** 2)
    ) ** 0.5

Clad=1


core = box(0, 0, 0.5, 0.39) #500x390nm
polygons = OrderedDict(
    core=core,
    box=clip_by_rect(core.buffer(1.5, resolution=4), -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(core.buffer(1.5, resolution=4), -np.inf, 0, np.inf, np.inf),
)

resolutions = {"core": {"resolution": 0.025, "distance": 2.}}

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.6))

num_modes = 2

#For SFWM we have energy conservation and momemtum(k) conservation for 2pumps and signal+idler
lambda_p0 = 1.4
lambda_s0 = 1.097
lambda_i0 = 1.686

basis0 = Basis(mesh, ElementTriP0())

epsilon_p = basis0.zeros()
epsilon_s = basis0.zeros()
epsilon_i = basis0.zeros()


for wavelength, epsilon in zip([lambda_p0, lambda_s0, lambda_i0], [epsilon_p, epsilon_s, epsilon_i]):
    for subdomain, n_func in {
        "core": n_X,
        "box": n_silicon_dioxide,
        "clad": lambda _: Clad,
    }.items():
        n = n_func(wavelength)
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2


modes_p = compute_modes(basis0, epsilon_p, wavelength=lambda_p0, num_modes=num_modes, order=1)
modes_s = compute_modes(basis0, epsilon_s, wavelength=lambda_s0, num_modes=num_modes, order=1)
modes_i = compute_modes(basis0, epsilon_i, wavelength=lambda_i0, num_modes=num_modes, order=1)


mode_p = max(modes_p, key=lambda mode: mode.te_fraction)
mode_s = max(modes_s, key=lambda mode: mode.te_fraction)
mode_i = max(modes_i, key=lambda mode: mode.te_fraction)

A_eff = calculate_sfwm_Aeff(basis0, mode_p, mode_s, mode_i)
print(A_eff)

#Calculation for non-linear coef
chi_3 = 5e-21  # m^2/V^2  #7e-20?
lambda_p0_m = lambda_p0 * 1e-6  #to m
n_p0 = np.real(mode_p.n_eff)
A_eff_m2 = A_eff * 1e-12  #to m^2

omega_p0 = 2 * np.pi * c / lambda_p0_m

gamma = (3 * chi_3 * omega_p0) / (4 * epsilon_0 * c**2 * n_p0**2 * A_eff_m2)

print("gamma:",gamma)
