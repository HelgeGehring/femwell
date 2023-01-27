import matplotlib.pyplot as plt
import numpy as np
from skfem import (
    Basis,
    BilinearForm,
    ElementTriN2,
    ElementTriP0,
    ElementTriP2,
    FacetBasis,
    Mesh,
    asm,
    condense,
    solve,
)
from skfem.helpers import curl, dot, grad, inner
from tqdm import tqdm

from femwell.mode_solver import solver_slepc


def compute_modes(basis_epsilon_r, epsilon_r, mu_r, num_modes, phase_x):
    one_over_u_r = 1

    basis = basis_epsilon_r.with_element(ElementTriN2() * ElementTriP2())

    @BilinearForm(dtype=complex)
    def aform(E, lam, v, mu, w):
        return (
            one_over_u_r * curl(E) * curl(v)
            + one_over_u_r * dot(grad(lam), v)
            + w["epsilon"] * dot(E, grad(mu))
            - w["epsilon"] * lam * mu
        )

    @BilinearForm(dtype=complex)
    def bform(E, lam, v, mu, w):
        return w["epsilon"] * dot(E, v)

    A = aform.assemble(basis, epsilon=basis0.interpolate(epsilon_r))
    B = bform.assemble(basis, epsilon=basis0.interpolate(epsilon_r))

    @BilinearForm(dtype=complex)
    def penalty(u, u_, v, v_, w):
        u1 = (w.idx[0] == 0) * u
        u2 = (w.idx[0] == 1) * u
        v1 = (w.idx[1] == 0) * v
        v2 = (w.idx[1] == 1) * v
        ju = u1 - w["phase"] * u2
        jv = v1 - w["phase"] * v2

        u1 = (w.idx[0] == 0) * u_
        u2 = (w.idx[0] == 1) * u_
        v1 = (w.idx[1] == 0) * v_
        v2 = (w.idx[1] == 1) * v_
        ju_ = u1 - w["phase"] * u2
        jv_ = v1 - w["phase"] * v2

        return 1.0 / 1e-2 * inner(ju, jv) + 1.0 / 1e-2 * inner(ju_, jv_)

    fbases = [
        FacetBasis(basis.mesh, basis.elem, facets="left"),
        FacetBasis(basis.mesh, basis.elem, facets="right"),
    ]
    D1 = asm(penalty, fbases, fbases, phase=phase_x)

    fbases = [
        FacetBasis(basis.mesh, basis.elem, facets="top"),
        FacetBasis(basis.mesh, basis.elem, facets="bottom"),
    ]
    D2 = asm(penalty, fbases, fbases, phase=1)

    lams, xs = solve(
        *condense(A + D1, B, D=basis.get_dofs(["top", "bottom"]), x=basis.zeros(dtype=complex)),
        solver=solver_slepc(k=num_modes, which="LR", sigma=10)
    )

    return np.sqrt(lams), basis, xs


if __name__ == "__main__":
    from collections import OrderedDict

    from mesh import mesh_from_OrderedDict
    from shapely.geometry import LineString, box

    width = 0.4
    cell_width = 0.45
    cell_height = 4

    structure = box(0, -width / 2, cell_width, width / 2)
    hole = structure.centroid.buffer(0.1)
    cell = box(0, -cell_height / 2, cell_width, cell_height / 2) - structure

    polygons = OrderedDict(
        left=LineString(((0, -cell_height / 2), (0, cell_height / 2))),
        right=LineString(((cell_width, -cell_height / 2), (cell_width, cell_height / 2))),
        top=LineString(((0, cell_height / 2), (cell_width, cell_height / 2))),
        bottom=LineString(((0, -cell_height / 2), (cell_width, -cell_height / 2))),
        hole=hole,
        structure=structure,
        cell=cell.geoms[0],
        cell1=cell.geoms[1],
    )

    resolutions = dict(
        hole={"resolution": 0.05, "distance": 1},
        core={"resolution": 0.05, "distance": 1},
    )

    mesh = mesh_from_OrderedDict(
        polygons, resolutions, filename="mesh.msh", default_resolution_max=0.05
    )
    mesh = Mesh.load("mesh.msh")

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros(dtype=complex) + 1**2
    epsilon[basis0.get_dofs(elements="structure")] = 3**2
    # basis0.plot(np.real(epsilon), colorbar=True).show()

    phases = np.linspace(0, 0.5, 5) * cell_width
    # phases = [.4*cell_width]
    print(phases)
    results = []
    for phase in tqdm(phases):
        lams, basis, xs = compute_modes(
            basis0, epsilon, 1, 5, phase_x=np.exp(2j * np.pi * phase / cell_width)
        )
        results.append((lams[:5], basis, xs))

    lams = np.array([np.sort(result[0]) / (2 * np.pi / cell_width) for result in results])
    plt.plot(phases, lams)
    plt.show()

    from mode_solver import plot_mode

    plot_mode(basis, np.real(xs[:, -3]), direction="x", colorbar="same")
    plt.show()

    plot_mode(basis, np.imag(xs[:, -3]), direction="x", colorbar="same")
    plt.show()
