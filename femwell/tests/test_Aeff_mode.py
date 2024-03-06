from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
from shapely.ops import clip_by_rect
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from shapely.geometry import box
from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict


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

core = box(0, 0, 0.5, 0.39)
polygons = OrderedDict(
    core=core,
    box=clip_by_rect(core.buffer(1.5, resolution=4), -np.inf, -np.inf, np.inf, 0),
    clad=clip_by_rect(core.buffer(1.5, resolution=4), -np.inf, 0, np.inf, np.inf),
)

resolutions = {"core": {"resolution": 0.025, "distance": 2.}}

mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.6))

num_modes = 2

lambda_p0 = 1.4
lambda_s0 = 1.097
lambda_i0 = 1.686

basis0 = Basis(mesh, ElementTriP0())

epsilon_p = basis0.zeros(dtype=complex)
epsilon_s = basis0.zeros(dtype=complex)
epsilon_i = basis0.zeros(dtype=complex)

for wavelength, epsilon in zip([lambda_p0, lambda_s0, lambda_i0], [epsilon_p, epsilon_s, epsilon_i]):
    for subdomain, n_func in {
        "core": n_X,
        "box": n_silicon_dioxide,
        "clad": lambda _: Clad,
    }.items():
        n = n_func(wavelength)
        epsilon[basis0.get_dofs(elements=subdomain)] = n**2
modes = compute_modes(basis0, epsilon, wavelength=wavelength, num_modes=2, order=2)


modes_p = compute_modes(basis0, epsilon_p, wavelength=lambda_p0, num_modes=num_modes, order=1)
modes_s = compute_modes(basis0, epsilon_s, wavelength=lambda_s0, num_modes=num_modes, order=1)
modes_i = compute_modes(basis0, epsilon_i, wavelength=lambda_i0, num_modes=num_modes, order=1)

mode_p = max(modes_p, key=lambda mode: mode.te_fraction)
mode_s = max(modes_s, key=lambda mode: mode.te_fraction)
mode_i = max(modes_i, key=lambda mode: mode.te_fraction)

print(f"Effective refractive index for p mode: {np.real(mode_p.n_eff):.4f}")
mode_p.show(mode_p.E.real, colorbar=True, direction="x")
mode_p.show(mode_p.E.imag, colorbar=True, direction="x")
fig, ax = plt.subplots()
mode_p.plot_intensity(ax=ax)
plt.title("Normalized Intensity for p mode")
plt.tight_layout()
plt.show()

print(f"Effective refractive index for s mode: {np.real(mode_s.n_eff):.4f}")
mode_s.show(mode_s.E.real, colorbar=True, direction="x")
mode_s.show(mode_s.E.imag, colorbar=True, direction="x")
fig, ax = plt.subplots()
mode_s.plot_intensity(ax=ax)
plt.title("Normalized Intensity for s mode")
plt.tight_layout()
plt.show()

print(f"Effective refractive index for i mode: {np.real(mode_i.n_eff):.4f}")
mode_i.show(mode_i.E.real, colorbar=True, direction="x")
mode_i.show(mode_i.E.imag, colorbar=True, direction="x")
fig, ax = plt.subplots()
mode_i.plot_intensity(ax=ax)
plt.title("Normalized Intensity for i mode")
plt.tight_layout()
plt.show()
