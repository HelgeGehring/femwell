import matplotlib.pyplot as plt
import numpy as np
import shapely
from meshwell.model import Model
from meshwell.polysurface import PolySurface
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
    enforce,
    solve,
)
from skfem.helpers import curl, dot, grad, inner
from skfem.io import from_meshio
from skfem.io.meshio import from_meshio
from skfem.models.general import curluv
from tqdm import tqdm

from femwell.coulomb import solve_coulomb
from femwell.mesh import mesh_from_OrderedDict
from femwell.visualization import plot_domains, plot_subdomain_boundaries


def solve_magnetostatic_2D(basis_mu_r, mu_r, mesh, current_densities):
    """The 2D magnetostatic problem can be reduced to a scalar Coulomb equation with source."""

    basis = basis_mu_r.with_element(ElementTriP1())

    @BilinearForm
    def dudv(A, v, w):
        return w["mu_r"] * inner(grad(A), grad(v))

    @LinearForm
    def unit_load(v, _):
        return v

    for (
        i,
        (domain, current_density),
    ) in enumerate(current_densities.items()):
        basis_source = FacetBasis(mesh, basis.elem, facets=mesh.boundaries[domain], intorder=4)
        basis_source_1d = basis_source.with_element(ElementTriP1())
        if i == 0:
            J_rhs = basis_source_1d.zeros()
        J_rhs += current_density * unit_load.assemble(basis_source_1d)

    A = dudv.assemble(basis, mu_r=basis_mu_r.interpolate(mu_r))
    A, J_rhs = enforce(A, J_rhs, D=mesh.boundary_nodes())

    return basis, solve(A, J_rhs)


if __name__ == "__main__":

    # Define some parameters for the solenoid:
    core_permeability = 1
    diameter = 2
    length = 1
    current = 1

    def continuous_solenoid_cross_section_mesh(
        diameter=diameter,
        length=length,
        dx_domain=10,
        dy_domain=10,
    ):

        solenoid_polygon = shapely.box(-diameter / 2, -length / 2, diameter / 2, length / 2)
        left_air_polygon = shapely.box(-dx_domain / 2, -dy_domain / 2, -diameter / 2, dy_domain / 2)
        right_air_polygon = shapely.box(diameter / 2, -dy_domain / 2, dx_domain / 2, dy_domain / 2)
        center_air_polygon = shapely.box(-diameter / 2, -dy_domain / 2, diameter / 2, dy_domain / 2)

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
                verbosity=0,
            )
        )

    mesh = continuous_solenoid_cross_section_mesh(
        diameter=diameter,
        length=length,
        dx_domain=10,
        dy_domain=10,
    )
    # mesh.draw().show()

    # plot_domains(mesh)
    # plt.show()

    # plot_subdomain_boundaries(mesh)

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
                "solenoid___left": current,
                "solenoid___right": -current,
            },
        )

        return basis_A, A, basis_mu_r, mu_r

    basis_A, A, basis_mu_r, mu_r = magnetic_potential(
        mesh, current=current, permeability_mu_r=core_permeability
    )

    fig, ax = plt.subplots()
    for subdomain in ["solenoid"]:
        basis_mu_r.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
    basis_A.plot(A, ax=ax, shading="gouraud", colorbar=True)

    bfield = asm(
        LinearForm(curluv).partial(basis_A.interpolate(A)),
        basis_A.with_element(ElementVector(ElementTriP1())),
    )

    print(np.min(basis_A.doflocs[0, :]))
    print(np.max(basis_A.doflocs[0, :]))
    print(np.min(basis_A.doflocs[1, :]))
    print(np.max(basis_A.doflocs[1, :]))

    # Define the grid for interpolation
    Nx, Ny = 50, 50
    xmin, xmax = mesh.p[0].min(), mesh.p[0].max()
    ymin, ymax = mesh.p[1].min(), mesh.p[1].max()
    grid_x, grid_y = np.meshgrid(np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny))

    # Interpolate the magnetic field on the grid
    from scipy.interpolate import griddata

    bfield_x = bfield.reshape((-1, 2))[:, 0]
    bfield_y = bfield.reshape((-1, 2))[:, 1]
    grid_bfield_x = griddata(mesh.p.T, bfield_x, (grid_x, grid_y), method="cubic")
    grid_bfield_y = griddata(mesh.p.T, bfield_y, (grid_x, grid_y), method="cubic")

    # Plot the streamlines
    ax.streamplot(grid_x, grid_y, grid_bfield_x, grid_bfield_y, color="black", linewidth=0.5)

    plt.show()

    # curl_A = curl(basis_A.interpolate(A))
    # print(A.shape)
    # print(curl_A.shape)
    # curl_A_x = basis_A.interpolate(curl_A[0,:])
    # curl_A_y = basis_A.interpolate(curl_A[1,:])
    # coordinates = basis_A.doflocs

    # velocity = asm(
    #     LinearForm(curluv).partial(ib.interpolate(psi)),
    #     ib.with_element(ElementVector(ElementTriP1())),
    # )

    # fig, ax = plt.subplots()
    # ax.quiver(coordinates[0, :], coordinates[1, :], curl_A[0,:,0], curl_A[1,:,0])
    # plt.show()
