# implementing https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-6-1053&id=312806

from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

from skfem import *
from skfem.helpers import *
from skfem.io import from_meshio
import shapely

from mesh import mesh_from_OrderedDict
from solver import solver_eigen_slepc

height = 1
a = .330
b = .7
c = .2

wavelength = .8
k0 = 2*np.pi/wavelength
k0 = .7/a
print(k0)

print(k0, k0*a)

left = shapely.LineString([(0, y) for y in np.linspace(-height, height, 2)])
right = shapely.LineString([(a, y) for y in np.linspace(-height, height, 2)])
top = shapely.LineString([(x,height) for x in np.linspace(0, a, 2)])
bottom = shapely.LineString([(x,-height) for x in np.linspace(0, a, 2)])

box = shapely.box(0,-height,a,height)
structure = shapely.box(0,-b/2,a,b/2)
hole = shapely.box(a/4,-c/2,a/4*3,c/2)

resolutions = {
    'box': {'resolution':.1, 'distance':1}
}

mesh = from_meshio(mesh_from_OrderedDict(OrderedDict(
    left=left, right=right, top=top, bottom=bottom,
    hole=hole,
    structure=structure,
    box=box,
), resolutions=resolutions, filename='mesh.msh', default_resolution_max=.05))


basis_vec = Basis(mesh, ElementTriP1()*ElementTriP1())
basis0 = basis_vec.with_element(ElementTriP0())
basis1 = basis_vec.with_element(ElementTriP1())

epsilon = basis0.zeros(dtype=np.complex64) + 1.45
epsilon[basis0.get_dofs(elements='structure')] = 3.5
epsilon**=2

# basis0.plot(np.real(epsilon), ax=basis1.draw(),colorbar=True).show()

@BilinearForm(dtype=np.complex64)
def A(phi, k_phi, v, k_v, w):
    return -dot(grad(phi), grad(v)) + k0**2 * w.epsilon * phi * v

@BilinearForm(dtype=np.complex64)
def B(phi, k_phi, v, k_v, w):
    return 2j * grad(k_phi)[0] * v

@BilinearForm(dtype=np.complex64)
def C(phi, k_phi, v, k_v, w):
    return - k_phi * v

@BilinearForm(dtype=np.complex64)
def I_k_phi(phi, k_phi, v, k_v, w):
    return k_phi * k_v

@BilinearForm(dtype=np.complex64)
def I_phi(phi, k_phi, v, k_v, w):
    return phi * k_v

fbases = [
    FacetBasis(mesh, basis_vec.elem, facets='left'),
    FacetBasis(mesh, basis_vec.elem, facets='right'),
]

assert np.all(fbases[0].default_parameters()['x'][1] == fbases[1].default_parameters()['x'][1])


@BilinearForm(dtype=np.complex64)
def penalty(phi, k_phi, v, k_v, w):
    u1 = (w.idx[0] == 0) * phi
    u2 = (w.idx[0] == 1) * phi
    v1 = (w.idx[1] == 0) * v
    v2 = (w.idx[1] == 1) * v
    ju = u1 - u2
    jv = v1 - v2

    u1 = (w.idx[0] == 0) * k_phi
    u2 = (w.idx[0] == 1) * k_phi
    v1 = (w.idx[1] == 0) * k_v
    v2 = (w.idx[1] == 1) * k_v
    k_ju = u1 - u2
    k_jv = v1 - v2

    return 1. / 1e-2 * (ju * jv + k_ju * k_jv)

A = A.assemble(basis_vec, epsilon=basis0.interpolate(epsilon)) + B.assemble(basis_vec) + I_k_phi.assemble(basis_vec) + asm(penalty, fbases, fbases)
f = - C.assemble(basis_vec) + I_phi.assemble(basis_vec)

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

def solver_dense(**kwargs):
    def solver(A,B):
        import scipy.linalg
        return scipy.linalg.eig(A.todense(),B.todense())
    return solver

ks, xs = solve(
    *condense(
        A,
        f,
        D=basis_vec.get_dofs({'top','bottom'}),
        #x=basis_vec.zeros(dtype=np.complex64),
        expand=True
    ),
    #solver=solver_eigen_slepc(k=15, which='SM', sigma=k0/10)
    solver=solver_dense(k=5, which='LM', sigma=np.pi/2)
)

xs = xs[:basis_vec.N]

#idx = np.abs(ks) < 1e10
#ks = ks[idx]
#xs = xs[:,idx]

idx = np.abs(np.abs(ks)).argsort()[::1]   
ks = ks[idx]
xs = xs[:,idx]

print(ks)

plt.plot(np.real(ks))
plt.show()

(phis, basis_phi), (k_phis, basis_k_phi) = basis_vec.split(xs)

for i in range(phis.shape[-1]):
    ax=basis1.draw()
    ax.set_title(f'{ks[i]}')
    basis_phi.plot(np.real(phis[...,i]), ax=ax, shading='gouraud', colorbar=True).show()
basis_phi.plot(np.imag(phis[...,150]), ax=basis1.draw(), shading='gouraud').show()