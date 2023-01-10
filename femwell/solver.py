import numpy as np

def solver_eigen_slepc(**kwargs):
    params = {
        'sigma': None,
        'k': 5,
        'which': 'TM',
        'st': None
    }

    params.update(kwargs)

    def solver(K, M, **solve_time_kwargs):
        from petsc4py import PETSc
        from slepc4py import SLEPc
        which = {
            'LM': SLEPc.EPS.Which.LARGEST_MAGNITUDE,
            'SM': SLEPc.EPS.Which.SMALLEST_MAGNITUDE,
            'LR': SLEPc.EPS.Which.LARGEST_REAL,
            'SR': SLEPc.EPS.Which.SMALLEST_REAL,
            'LI': SLEPc.EPS.Which.LARGEST_IMAGINARY,
            'SI': SLEPc.EPS.Which.SMALLEST_IMAGINARY,
            'TM': SLEPc.EPS.Which.TARGET_MAGNITUDE,
            'TR': SLEPc.EPS.Which.TARGET_REAL,
            'TI': SLEPc.EPS.Which.TARGET_IMAGINARY,
            'ALL': SLEPc.EPS.Which.ALL,
        }
        st = {
            'CAYLEY': SLEPc.ST.Type.CAYLEY,
            'FILTER': SLEPc.ST.Type.FILTER,
            'PRECOND': SLEPc.ST.Type.PRECOND,
            'SHELL': SLEPc.ST.Type.SHELL,
            'SHIFT': SLEPc.ST.Type.SHIFT,
            'SINVERT': SLEPc.ST.Type.SINVERT
        }

        params.update(solve_time_kwargs)


        filled_spots = set(zip(K.tocoo().row, K.tocoo().col))
        for i in range(K.shape[0]):
            if (i,i) not in filled_spots:
                K[i,i]=0
        filled_spots = set(zip(M.tocoo().row, M.tocoo().col))
        for i in range(K.shape[0]):
            if (i,i) not in filled_spots:
                M[i,i]=0

        K_ = PETSc.Mat().createAIJ(size=K.shape, csr=(K.indptr, K.indices, K.data))
        M_ = PETSc.Mat().createAIJ(size=M.shape, csr=(M.indptr, M.indices, M.data))


        eps = SLEPc.EPS().create()
        eps.setDimensions(params['k'])
        eps.setOperators(K_, M_)
        eps.setType(SLEPc.EPS.Type.LAPACK)
        if params['st']:
            eps.getST().setType(st[params['st']])
        eps.setWhichEigenpairs(which[params['which']])
        if params['sigma']:
            eps.setTarget(params['sigma'])
        eps.solve()

        xr, xi = K_.getVecs()
        lams, xs = [], []
        for i in range(eps.getConverged()):
            val = eps.getEigenpair(i, xr, xi)
            lams.append(val)
            xs.append(np.array(xr) + 1j * np.array(xi))

        return np.array(lams), np.array(xs).T

    return solver

from petsc4py import PETSc
from slepc4py import SLEPc
print(dir(SLEPc.EPS.Type))