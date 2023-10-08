# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: waveguidesmodes
#     language: python
#     name: python3
# ---

# # Parallel plate capacitor
#
# The capacitance $C$ of a system given a potential difference $\Delta V$ between two conductors can be calculated from
#
# $$ C = \frac{2W}{(\Delta V)^2} $$
#
# with
#
# $$ W = \frac{1}{2} \int_\Omega \epsilon E \cdot E d\Omega $$
#
# where $\epsilon$ is the permittivity distribution, $E$ the electric field, and $\Omega$ the domain. The integrand is only non-zero close to the conductors, where the field is concentrated.
#
# In this notebook we compute the capacitance of a parallel-plate capacitor, and show that we recover the theoretical result $C = \frac{\epsilon A}{d}$ in the limit of negligible fringe fields.
#
# First, we parametrize a simple geometry:

# +
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.affinity
from scipy.constants import epsilon_0, speed_of_light
from shapely.ops import clip_by_rect
from skfem import Basis, ElementDG, ElementTriP0, ElementTriP1
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict
from femwell.visualization import plot_domains, plot_subdomain_boundaries

from meshwell.model import Model
from meshwell.polysurface import PolySurface

from femwell.coulomb import solve_coulomb

from skfem.helpers import dot
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
# -

dielectric_epsilon = 16
separation = 4


def parallel_plate_capacitor_mesh(width, 
                                  separation=separation,
                                    thickness=1,
                                    ):

    top_plate_polygon = shapely.geometry.box(-width / 2, separation / 2, width / 2, separation / 2 + thickness)
    bottom_plate_polygon = shapely.geometry.box(-width / 2, -separation / 2 - thickness, width / 2, -separation / 2)
    dielectric_polygon = shapely.geometry.box(-width / 2, -separation / 2, width / 2, separation / 2)

    capacitor_polygon = shapely.unary_union([top_plate_polygon, bottom_plate_polygon, dielectric_polygon])
    air_polygon = capacitor_polygon.buffer(20, resolution=8)


    model = Model()

    top_plate = PolySurface(
        polygons=top_plate_polygon, 
        model=model, 
        physical_name="top_plate", 
        mesh_bool=False,
        resolution={"resolution": 0.5, "DistMax": 2},
        mesh_order=1
    )
    bottom_plate = PolySurface(
        polygons=bottom_plate_polygon, 
        model=model, 
        physical_name="bottom_plate", 
        mesh_bool=False,
        resolution={"resolution": 0.5, "DistMax": 2},
        mesh_order=2
    )
    dielectric = PolySurface(
        polygons=dielectric_polygon, 
        model=model, 
        physical_name="dielectric", 
        mesh_bool=True,
        mesh_order=3
    )
    air = PolySurface(
        polygons=air_polygon, 
        model=model, 
        physical_name="air", 
        mesh_bool=True,
        mesh_order=4
    )

    return from_meshio(model.mesh(entities_list=[top_plate, bottom_plate, dielectric, air], 
                      filename="mesh.msh",
                      default_characteristic_length=0.5,
                      ))

mesh = parallel_plate_capacitor_mesh(width = 2, 
                                  separation = 2,
                                    thickness=1,
                                    )
mesh.draw().show()

plot_domains(mesh)
plt.show()

plot_subdomain_boundaries(mesh)


# Since the electric field inside a metal is zero, we define the voltage on the boundary of the metal, and exclude the metal from the mesh:

# +
def potential(mesh, dV = 1, dielectric_epsilon = 16):

    basis_epsilon = Basis(mesh, ElementTriP0())
    epsilon = basis_epsilon.ones()

    epsilon[basis_epsilon.get_dofs(elements=("dielectric"))] = dielectric_epsilon

    basis_u, u = solve_coulomb(
        basis_epsilon,
        epsilon,
        {"top_plate___dielectric": dV, "top_plate___air": dV, "bottom_plate___dielectric": 0, "bottom_plate___air": 0,}
    )

    return basis_u, u, basis_epsilon, epsilon

basis_u, u, basis_epsilon, epsilon = potential(mesh, dV=1)
# -

fig, ax = plt.subplots()
for subdomain in basis_epsilon.mesh.subdomains.keys() - {"gmsh:bounding_entities"}:
    basis_epsilon.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
basis_u.plot(u, ax=ax, shading="gouraud", colorbar=True)
# basis_vec.plot(-u_grad, ax=ax)
plt.show()

# +
for subdomain in basis_epsilon.mesh.subdomains.keys() - {"gmsh:bounding_entities"}:
    basis_epsilon.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
basis_grad = basis_u.with_element(ElementDG(basis_u.elem))

fig, ax = plt.subplots()
e_x = basis_u.project(-basis_epsilon.interpolate(epsilon) * basis_u.interpolate(u).grad[0])
basis_u.plot(e_x, ax=ax, shading="gouraud", colorbar=True)
plt.show()

fig, ax = plt.subplots()
e_y = basis_u.project(-basis_epsilon.interpolate(epsilon) * basis_u.interpolate(u).grad[1])
basis_u.plot(e_y, ax=ax, shading="gouraud", colorbar=True)
plt.show()


# -

def capacitance(width = 2, 
                separation = separation,
                    thickness=1,
                    dV = 2,
                    dielectric_epsilon = dielectric_epsilon,
                    ):

    mesh = parallel_plate_capacitor_mesh(width = width, 
                                    separation = separation,
                                        thickness=thickness,
                                        )
    basis_u, u, basis_epsilon, epsilon = potential(mesh, dV=dV, dielectric_epsilon = dielectric_epsilon)
    
    @Functional(dtype=complex)
    def W(w):
        return 0.5 * w["epsilon"] * dot(w["u"].grad, w["u"].grad)

    C = 2 * W.assemble(
        basis_u,
        epsilon=basis_epsilon.interpolate(epsilon),
        u=basis_u.interpolate(u),
    ) / dV**2

    return C


widths = np.linspace(1, 50, 21)
Cs = []
for width in widths:
    Cs.append(capacitance(width = width, 
                separation = 4,
                    thickness=2,
                    dV = 2,
                    ))

plt.plot(widths, np.array(Cs))
plt.plot(widths, np.array(widths) * dielectric_epsilon / separation)


