# implementing https://opg.optica.org/ol/fulltext.cfm?uri=ol-40-6-1053&id=312806

from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

from skfem import *
from skfem.helpers import *
from skfem.utils import mpc
from skfem.io import from_meshio
import shapely

from mesh import mesh_from_OrderedDict
from solver import solver_eigen_scipy_operator

height = 5.76/2+5
a = .100
b = .78
c = .2
slab = .920+5
pml = 3

wavelength = 1
k0 = 2*np.pi/wavelength

left = shapely.LineString([(0, y) for y in np.linspace(-height, height, 2)])
right = shapely.LineString([(a, y) for y in np.linspace(-height, height, 2)])
top = shapely.LineString([(x,height) for x in np.linspace(0, a, 2)])
bottom = shapely.LineString([(x,-height) for x in np.linspace(0, a, 2)])

box = shapely.box(0,-height,a,height)
structure = shapely.box(0,-b/2,a,b/2)
structure1 = shapely.box(0,height-slab,a,height)
structure2 = shapely.box(0,-height+slab,a,-height)

resolutions = {
    'structure': {'resolution':.1, 'distance':.1},
    'hole': {'resolution':.1, 'distance':.1}
}

mesh = from_meshio(mesh_from_OrderedDict(OrderedDict(
    left=left, right=right, top=top, bottom=bottom,
    structure=structure,
    structure1=structure1,
    structure2=structure2,
    box=box,
), resolutions=resolutions, filename='mesh.msh', default_resolution_max=.1))

basis_vec = Basis(mesh, ElementTriP1()*ElementTriP1())
basis0 = basis_vec.with_element(ElementDG(ElementTriP1()))

epsilon_r = basis0.zeros(dtype=np.complex64) + 1.45
epsilon_r[basis0.get_dofs(elements='box')] = 1.39
epsilon_r**=2
basis0.plot(np.real(epsilon_r), ax=mesh.draw(),colorbar=True).show()
basis0.plot(np.imag(epsilon_r), ax=mesh.draw(),colorbar=True).show()

epsilon_r += basis0.project(lambda x: (.2j)*(np.clip(np.abs(x[1])-height+pml, 0, np.inf)/pml)**2, dtype=np.complex64)

@BilinearForm(dtype=np.complex64)
def A(phi, k_phi, v, k_v, w):
    return -d(phi)[0] * d(v)[0] - d(phi)[1] * d(v)[1] + k0**2 * (w.epsilon) * phi * v

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

A = A.assemble(basis_vec, epsilon=basis0.interpolate(epsilon_r)) + B.assemble(basis_vec) + I_k_phi.assemble(basis_vec)
f = - C.assemble(basis_vec) + I_phi.assemble(basis_vec)

left = basis_vec.get_dofs(facets='left')
right = basis_vec.get_dofs(facets='right')
top = basis_vec.get_dofs(facets='top')
bottom = basis_vec.get_dofs(facets='bottom')

ks, xs = solve(*mpc(-A, -f, M=left, S=np.concatenate((right, top, bottom))), solver=solver_eigen_scipy_operator(k=10, which='LM', sigma=k0 * epsilon_r.real.max()))

print(ks)

plt.plot(np.real(ks))
plt.plot(np.imag(ks))
plt.show()

(phis, basis_phi), (k_phis, basis_k_phi) = basis_vec.split(xs)


for i in range(xs.shape[-1]):
    fig, ax = plt.subplots(1,1,figsize=(10,4))
    plt.title(f'{ks[i]}')
    vminmax = np.max(np.abs(basis_phi.interpolate(phis[...,i])))


    for i_plot in range(10):
        phases = basis_phi.project(lambda x: np.exp(1j*ks[i]*(x[0]+i_plot*a)), dtype=np.complex64)
        phi_with_phase = basis_phi.project(basis_phi.interpolate(phis[...,i])*basis_phi.interpolate(phases), dtype=np.complex64) 
        im = ax.tripcolor(mesh.p[0]+i_plot*a, mesh.p[1], mesh.t.T, np.real(phi_with_phase), cmap='seismic', shading='gouraud', vmin=-vminmax, vmax=vminmax)
    plt.show()