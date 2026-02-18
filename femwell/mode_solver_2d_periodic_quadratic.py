# implementing https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-6-1053&id=312806

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import shapely
from mesh import mesh_from_OrderedDict
from skfem import (
    Basis,
    BilinearForm,
    ElementTriP0,
    ElementTriP1,
    ElementTriP2,
    FacetBasis,
    asm,
)
from skfem.helpers import d, grad
from skfem.io import from_meshio

height = 5.76 / 2 + 5
a = 0.100
b = 0.78
c = 0.2
slab = 0.920 + 5
pml = 3

wavelength = 1
k0 = 2 * np.pi / wavelength
# k0 = .75/a
print(k0)

print(k0, k0 * a)

left = shapely.LineString([(0, y) for y in np.linspace(-height, height, 2)])
right = shapely.LineString([(a, y) for y in np.linspace(-height, height, 2)])
top = shapely.LineString([(x, height) for x in np.linspace(0, a, 2)])
bottom = shapely.LineString([(x, -height) for x in np.linspace(0, a, 2)])

box = shapely.box(0, -height, a, height)
structure = shapely.box(0, -b / 2, a, b / 2)
structure1 = shapely.box(0, height - slab, a, height)
structure2 = shapely.box(0, -height + slab, a, -height)
hole = shapely.box(a / 4, -c / 2, a / 4 * 3, c / 2)

resolutions = {
    "structure": {"resolution": 0.1, "distance": 0.1},
    "hole": {"resolution": 0.1, "distance": 0.1},
}

mesh = from_meshio(
    mesh_from_OrderedDict(
        OrderedDict(
            left=left,
            right=right,
            top=top,
            bottom=bottom,
            # hole=hole,
            structure=structure,
            structure1=structure1,
            structure2=structure2,
            box=box,
        ),
        resolutions=resolutions,
        filename="mesh.msh",
        default_resolution_max=0.1,
    )
)


basis_vec = Basis(mesh, ElementTriP1())
basis0 = basis_vec.with_element(ElementTriP0())
basis1 = basis_vec.with_element(ElementTriP1())
basis_pml = basis_vec.with_element(ElementTriP2())

epsilon = basis0.zeros(dtype=np.complex64) + 1.39
epsilon = basis0.project(lambda x: 1.45)
# epsilon[basis0.get_dofs(elements='structure')] = 1.45
# epsilon[basis0.get_dofs(elements='structure1')] = 1.45
# epsilon[basis0.get_dofs(elements='structure2')] = 1.45
epsilon[basis0.get_dofs(elements="box")] = 1.39
epsilon **= 2

pml = basis_pml.project(
    lambda x: 1 + 5j * (np.clip(np.abs(x[1]) - height + pml, 0, np.inf) / pml) ** 2,
    dtype=np.complex64,
)

basis0.plot(np.real(epsilon), ax=basis1.draw(), colorbar=True).show()
basis_pml.plot(np.imag(pml), ax=basis1.draw(), colorbar=True).show()


@BilinearForm(dtype=np.complex64)
def A(phi, v, w):
    return -d(phi)[0] * d(v)[0] - d(phi)[1] * d(v)[1] / w.pml**2 + k0**2 * w.epsilon * phi * v


@BilinearForm(dtype=np.complex64)
def B(phi, v, w):
    return 2j * grad(phi)[0] * v


@BilinearForm(dtype=np.complex64)
def C(phi, v, w):
    return -1 * phi * v


fbases = [
    FacetBasis(mesh, basis_vec.elem, facets="left"),
    FacetBasis(mesh, basis_vec.elem, facets="right"),
]

assert np.all(fbases[0].default_parameters()["x"][1] == fbases[1].default_parameters()["x"][1])


@BilinearForm(dtype=np.complex64)
def penalty(phi, v, w):
    u1 = (w.idx[0] == 0) * phi
    u2 = (w.idx[0] == 1) * phi
    v1 = (w.idx[1] == 0) * v
    v2 = (w.idx[1] == 1) * v
    ju = u1 - u2
    jv = v1 - v2

    return 1.0 / 1e-2 * (ju * jv)


fbases2 = [
    FacetBasis(mesh, basis_vec.elem, facets="top"),
    FacetBasis(mesh, basis_vec.elem, facets="bottom"),
]


@BilinearForm(dtype=np.complex64)
def penalty2(phi, v, w):
    u1 = (w.idx[0] == 0) * phi
    u2 = (w.idx[0] == 1) * phi
    v1 = (w.idx[1] == 0) * v
    v2 = (w.idx[1] == 1) * v
    ju = u1 - u2
    jv = v1 - v2

    return 1.0 / 1e-2 * (ju + jv)


from petsc4py import PETSc
from slepc4py import SLEPc

pep = SLEPc.PEP().create()

A = A.assemble(
    basis_vec, epsilon=basis0.interpolate(epsilon), pml=basis_pml.interpolate(pml)
) + asm(
    penalty, fbases, fbases
)  # +asm(penalty2, fbases2, fbases2)
B = B.assemble(basis_vec)
C = C.assemble(basis_vec)
mats = [PETSc.Mat().createAIJ(size=K.shape, csr=(K.indptr, K.indices, K.data)) for K in (A, B, C)]
pep.setOperators(mats)
pep.setDimensions(basis_vec.N * 2)

nev, ncv, mpd = pep.getDimensions()
print("")
print("Number of requested eigenvalues: %i" % nev)

pep.setTarget(1 / k0)
print("target", 1 / k0)
# pep.getST().setType((SLEPc.ST.Type.SINVERT))
pep.setWhichEigenpairs(SLEPc.PEP.Which.LARGEST_MAGNITUDE)
# pep.setType(SLEPc.PEP.Type.JD)
# pep.setProblemType(SLEPc.PEP.ProblemType.GENERAL)
print("set")


def monitor(eps, its, nconv, eig, err):
    print(its, nconv)


pep.setMonitor(monitor)
pep.solve()
print(nconv := pep.getConverged())

xr, xi = mats[0].createVecs()

ks = []
xs = []

for i in range(nconv):
    k = pep.getEigenpair(i, xr, xi)
    if np.abs(np.real(k)) < 8.8:
        continue
    error = pep.computeError(i)

    # print("%9f%+9f j    %12g" % (k.real, k.imag, error))
    ks.append(k)
    xs.append(np.array(xr))

ks = np.array(ks)
xs = np.array(xs).T
# print(np.array(xr))

# left = np.concatenate((basis_vec.get_dofs(facets='left').nodal['u^1'], basis_vec.get_dofs(facets='left').nodal['u^2']))
# right = np.concatenate((basis_vec.get_dofs(facets='right').nodal['u^1'], basis_vec.get_dofs(facets='right').nodal['u^2']))

# D = basis_vec.get_dofs({'top', 'bottom'}).all()

# left = np.setdiff1d(left, D)
# right = np.setdiff1d(right, D)

# B = np.zeros((len(left), basis_vec.N))
# B[range(len(left)), left] = 1
# B[range(len(right)), right] = -1

# A = bmat([[A, B.T],
#           [B, None]], 'csr')
# import scipy.sparse
# f = bmat([[f, scipy.sparse.csr_array(np.empty((f.shape[0], len(left)), dtype=np.complex64))],
#           [scipy.sparse.csr_array(np.empty((len(left),f.shape[1]), dtype=np.complex64)), None]], 'csr')

xs = xs[: basis_vec.N]

# idx = np.abs(ks) < 1e10
# ks = ks[idx]
# xs = xs[:,idx]

idx = np.abs(np.abs(ks)).argsort()[::1]
ks = ks[idx]
xs = xs[:, idx]

plt.plot(np.real(ks))
plt.plot(np.imag(ks))
plt.show()

# phis, basis_phi), (k_phis, basis_k_phi) = basis_vec.split(xs)

for i in range(xs.shape[-1]):
    fig, axs = plt.subplots(1, 10, figsize=(10, 4))
    plt.title(f"{ks[i]}")
    # ax = basis1.draw()

    vminmax = np.max(np.abs(basis_vec.interpolate(xs[..., i])))
    for j, ax in enumerate(axs):
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_axis_off()
        # ax.set_aspect(1)
        basis_vec.mesh.draw(ax=ax, boundaries=True, boundaries_only=True)

        phases = basis_vec.project(
            lambda x: np.exp(1j * ks[i] * (x[0] + j * a)), dtype=np.complex64
        )

        # for subdomain in basis_vec.mesh.subdomains.keys() - {'gmsh:bounding_entities'}:
        #    basis_vec.mesh.restrict(subdomain).draw(ax=ax, boundaries_only=True)
        basis_vec.plot(
            np.imag(xs[..., i] * phases),
            shading="gouraud",
            ax=ax,
            vmin=-vminmax,
            vmax=vminmax,
            cmap="seismic",
        )
    fig.subplots_adjust(wspace=0, hspace=0)
    plt.show()
