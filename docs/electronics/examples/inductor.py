# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: waveguidesmodes
#     language: python
#     name: python3
# ---

# # Magnetostatics and inductance
#
# We follow the theory outlined in the documentation to compute the magnetic field of various currents distributions.

# +
import matplotlib.pyplot as plt
import numpy as np
import shapely
from meshwell.model import Model
from meshwell.polysurface import PolySurface
from scipy.interpolate import griddata
from skfem import (
    Basis,
    BilinearForm,
    ElementDG,
    ElementTriN1,
    ElementTriP0,
    ElementTriP1,
    ElementVector,
    FacetBasis,
    Functional,
    LinearForm,
    asm,
    condense,
    solve,
)
from skfem.helpers import curl, dot, grad, inner
from skfem.io import from_meshio
from skfem.io.meshio import from_meshio
from skfem.models.general import curluv
from tqdm import tqdm

from femwell.magnetostatic import solve_magnetostatic_2D
from femwell.mesh import mesh_from_OrderedDict
from femwell.visualization import plot_domains, plot_subdomain_boundaries

# %load_ext autoreload
# %autoreload 2


# +
core_permeability = 1
diameter = 2
length = 5
current = 1

dr_domain = 50
dz_domain = 50
# -


# Make a mesh
#
# <div class="alert alert-success">
# Note: below we use meshwell instead of femwell's built-in backend. You can install meshwell with `pip install meshwell`
# </div>


# + tags=["remove-stderr", "hide-output"]
def continuous_solenoid_cross_section_mesh(
    diameter=diameter,
    length=length,
    dr_domain=10,
    dz_domain=10,
):

    solenoid_polygon = shapely.box(-diameter / 2, -length / 2, diameter / 2, length / 2)
    left_air_polygon = shapely.box(-dr_domain / 2, -dz_domain / 2, -diameter / 2, dz_domain / 2)
    right_air_polygon = shapely.box(diameter / 2, -dz_domain / 2, dr_domain / 2, dz_domain / 2)
    center_air_polygon = shapely.box(-diameter / 2, -dz_domain / 2, diameter / 2, dz_domain / 2)

    model = Model()

    solenoid = PolySurface(
        polygons=solenoid_polygon,
        model=model,
        physical_name="solenoid",
        mesh_bool=True,
        resolution={"resolution": 0.1, "DistMax": 2, "SizeMax": 0.2},
        mesh_order=1,
    )
    left_air = PolySurface(
        polygons=left_air_polygon,
        model=model,
        physical_name="left",
        mesh_bool=True,
        mesh_order=2,
    )
    right_air = PolySurface(
        polygons=right_air_polygon,
        model=model,
        physical_name="right",
        mesh_bool=True,
        mesh_order=2,
    )
    center_air = PolySurface(
        polygons=center_air_polygon,
        model=model,
        physical_name="center",
        mesh_bool=True,
        mesh_order=2,
    )

    return from_meshio(
        model.mesh(
            entities_list=[solenoid, left_air, center_air, right_air],
            filename="mesh.msh",
            default_characteristic_length=0.5,
        )
    )


mesh = continuous_solenoid_cross_section_mesh(
    diameter=diameter,
    length=length,
    dr_domain=dr_domain,
    dz_domain=dz_domain,
)
mesh.draw().show()

plot_domains(mesh)
plt.show()

plot_subdomain_boundaries(mesh)
# -


# We use these three background volumes to easily tag the left and right sides of the solenoid to define currents


# +
def magnetic_potential(mesh, current=1, permeability_mu_r=1):

    # Permeability defined constant-wise to handle sharp discontinuities
    basis_mu_r = Basis(mesh, ElementTriP0())
    mu_r = basis_mu_r.ones()
    mu_r[basis_mu_r.get_dofs(elements=("solenoid"))] = permeability_mu_r

    basis_A, A = solve_magnetostatic_2D(
        basis_mu_r=basis_mu_r,
        mu_r=mu_r,
        mesh=mesh,
        current_densities={
            "solenoid___left": -current / length / 2,  # convert current to current density
            "solenoid___right": current / length / 2,
        },
    )

    return basis_A, A, basis_mu_r, mu_r


basis_A, A, basis_mu_r, mu_r = magnetic_potential(
    mesh, current=current, permeability_mu_r=core_permeability
)
# -

# The magnetic field is $\vec{B} = \nabla \times \vec{A}$, or in terms of test functions,

# +
vector_basis = Basis(mesh, ElementVector(ElementTriP1()))

# Need a minus sign here?
bfield = -1 * asm(
    LinearForm(curluv).partial(basis_A.interpolate(A)),
    vector_basis,
)

# +
fig, ax = plt.subplots()
for subdomain in ["solenoid"]:
    basis_mu_r.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
basis_A.plot(A, ax=ax, shading="gouraud", colorbar=True)

# Define the grid for streamplot interpolation
Nx, Ny = 100, 100
xmin, xmax = mesh.p[0].min(), mesh.p[0].max()
ymin, ymax = mesh.p[1].min(), mesh.p[1].max()
grid_x, grid_y = np.meshgrid(np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny))

# Interpolate the magnetic field on the grid
bfield_x = bfield.reshape((-1, 2))[:, 0]
bfield_y = bfield.reshape((-1, 2))[:, 1]
grid_bfield_x = griddata(mesh.p.T, bfield_x, (grid_x, grid_y), method="cubic")
grid_bfield_y = griddata(mesh.p.T, bfield_y, (grid_x, grid_y), method="cubic")

# Plot the streamlines
lw = (
    5
    * np.sqrt(grid_bfield_x**2 + grid_bfield_y**2)
    / np.max(np.sqrt(grid_bfield_x**2 + grid_bfield_y**2))
)
ax.streamplot(grid_x, grid_y, grid_bfield_x, grid_bfield_y, linewidth=lw, density=1.2, color="k")
plt.show()
# -

max_A_femwell = np.max(A)

min_A_femwell = np.min(A)

max_B_femwell = np.max(np.sqrt(grid_bfield_x**2 + grid_bfield_y**2))

min_B_femwell = np.min(np.sqrt(grid_bfield_x**2 + grid_bfield_y**2))

# We can compare this to [the analytical solution](https://en.wikipedia.org/wiki/Solenoid#Finite_continuous_solenoid):

# +
from scipy.special import ellipe, ellipk, elliprf, elliprj


# See https://github.com/scipy/scipy/issues/4452
def ellipp(n, m):
    if m >= 1:
        raise ValueError("m must be < 1")
    y = 1 - m
    rf = elliprf(0, y, 1)
    rj = elliprj(0, y, 1, 1 - n)
    return rf + rj * n / 3


def ellippinc(phi, n, m):
    nc = np.floor(phi / np.pi + 0.5)
    phi -= nc * np.pi
    sin_phi = np.sin(phi)
    sin2_phi = sin_phi * sin_phi
    sin3_phi = sin2_phi * sin_phi
    x = 1 - sin2_phi
    y = 1 - m * sin2_phi
    rf = elliprf(x, y, 1)
    rj = elliprj(x, y, 1, 1 - n * sin2_phi)
    val = sin_phi * rf + sin3_phi * rj * n / 3
    if nc != 0:
        rp = ellipp(n, m)
        val += 2 * nc * rp
    return val


def m(rho, R, eta):
    return 4 * R * rho / ((R + rho) ** 2 + eta**2)


def n(rho, R):
    return 4 * R * rho / (R + rho) ** 2


def continuous_solenoid_analytical_A(rho, z, R, l):
    """For mu0 = 1, I = 1"""

    def term(eta):
        sqrt_term = np.sqrt((R + rho) ** 2 + eta**2)
        m_val = m(rho, R, eta)
        n_val = n(rho, R)
        k_term = ellipk(m_val)
        e_term = ellipe(m_val)
        inc_term = np.vectorize(ellippinc)(np.pi / 2, n_val, m_val)

        return (
            eta
            / sqrt_term
            * (
                (m_val + n_val - m_val * n_val) / (m_val * n_val) * k_term
                - 1 / m_val * e_term
                + (n_val - 1) / n_val * inc_term
            )
        )

    eta_plus = z + l / 2
    eta_minus = z - l / 2

    return R / (np.pi * l) * (term(eta_plus) - term(eta_minus))


"""
Conversions

R = a
h2 = n
k2 = m
"""


def continuous_solenoid_analytical_Bz(rho, z, R, l):
    """For mu0 = 1, I = 1"""

    def term(eta):
        m_val = m(rho, R, eta)
        n_val = n(rho, R)
        k_term = ellipk(m_val)
        inc_term = np.vectorize(ellippinc)(np.pi / 2, n_val, m_val)

        return eta * np.sqrt(m_val) * (k_term + (R - rho) / (R + rho) * inc_term)

    eta_plus = z + l / 2
    eta_minus = z - l / 2

    return 1 / (4 * np.pi * l * np.sqrt(R * rho)) * (term(eta_plus) - term(eta_minus))


def continuous_solenoid_analytical_Br(rho, z, R, l):
    """For mu0 = 1, I = 1"""

    def term(eta):
        m_val = m(rho, R, eta)
        k_term = ellipk(m_val)
        e_term = ellipe(m_val)

        return (m_val - 2) / np.sqrt(m_val) * k_term + 2 / np.sqrt(m_val) * e_term

    eta_plus = z + l / 2
    eta_minus = z - l / 2

    return np.sqrt(R) / (2 * np.pi * l * np.sqrt(rho)) * (term(eta_plus) - term(eta_minus))


zs = np.linspace(-dz_domain / 2, dz_domain / 2, 200)
rhos = np.linspace(-dr_domain / 2, dr_domain / 2, 200)
rhos_positive = rhos[rhos >= 0]

X, Y = np.meshgrid(rhos, zs)
X_positive, Y_positive = np.meshgrid(rhos_positive, zs)
A = continuous_solenoid_analytical_A(X, Y, R=diameter / 2, l=length)

Bz_positive = continuous_solenoid_analytical_Bz(X_positive, Y_positive, R=diameter / 2, l=length)
Br_positive = continuous_solenoid_analytical_Br(X_positive, Y_positive, R=diameter / 2, l=length)

fig, ax = plt.subplots()
stream = ax.streamplot(
    X_positive, Y_positive, Br_positive, Bz_positive, color="black", linewidth=0.5, density=0.75
)
ax.set_xlabel("rho")
ax.set_ylabel("z")
ax.set_title("Streamplot of Magnetic Field (Br, Bz) in Continuous Solenoid")

c = ax.contourf(X, Y, A, levels=50, cmap="jet")
fig.colorbar(c, ax=ax)
ax.set_xlabel("r")
ax.set_ylabel("z")
ax.set_title(r"Analytical $A_{\phi}$ of Continuous Solenoid")
plt.show()
# -

min_A_analytical = np.min(A)
max_A_analytical = np.max(A)
min_B_analytical = np.min(np.sqrt(Br_positive**2 + Bz_positive**2))
max_B_analytical = np.max(np.sqrt(Br_positive**2 + Bz_positive**2))

min_A_analytical, min_A_femwell

max_A_analytical, max_A_femwell

max_B_analytical, max_B_femwell

# # Inductance
#
# The inductance $L$ of a system given a current $I$ in a conductor can be calculated from
#
# $$ L = \frac{2W}{I^2} $$
#
# with
#
# $$ W = \frac{1}{2} \int_\Omega B \cdot H d\Omega = \frac{1}{2} \int_\Omega \mu H \cdot H d\Omega $$
#
# where $\mu$ is the permeability distribution, $B$ the magnetic field, $H = \mu^{-1} B$ the magnetization field, and $\Omega$ the domain. The integrand is only non-zero close to the conductors, where the field is concentrated.

# We can compare the analytical solution to the FEM solution:
