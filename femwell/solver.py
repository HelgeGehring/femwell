import numpy as np

def solver_dense(**kwargs):
    def solver(A, B):
        import scipy.linalg
        return scipy.linalg.eig(A.todense(), B.todense())
    return solver

def solver_eigen_scipy_operator(**kwargs):
    """Solve generalized eigenproblem using SciPy (ARPACK).

    Returns
    -------
    EigenSolver
        A solver function that can be passed to :func:`solve`.

    """
    params = {
        'sigma': 10,
        'k': 5,
    }
    params.update(kwargs)

    def solver(K, M, **solve_time_kwargs):
        params.update(solve_time_kwargs)
        from scipy.sparse.linalg import eigs
        from scipy.sparse.linalg import LinearOperator
        from scipy.sparse.linalg import factorized

        M_inv = factorized(M)
        ks, xs = eigs(LinearOperator(K.shape, matvec=lambda v: M_inv(K @ v), dtype=np.complex64) , **params)

        idx = np.abs(np.real(ks)).argsort()[::-1]   
        ks = ks[idx]
        xs = xs[:, idx]

        return ks, xs

    return solver

def solver_eigen_slepc(**kwargs):
    params = {
        'sigma': None,
        'k': 5,
        'which': 'TM',
        'st': 'SINVERT'
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


        K_ = PETSc.Mat().createAIJ(size=K.shape, csr=(K.indptr, K.indices, K.data))
        M_ = PETSc.Mat().createAIJ(size=M.shape, csr=(M.indptr, M.indices, M.data))


        eps = SLEPc.EPS().create()
        eps.setDimensions(params['k'])
        eps.setOperators(K_, M_)
        eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
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

if __name__ == '__main__':
    from petsc4py import PETSc
    from slepc4py import SLEPc

    import scipy.sparse

    pep = SLEPc.PEP().create()

    A = scipy.sparse.csr_array(([24.], ([0], [0])), shape=(1, 1), dtype=np.complex64)
    B= scipy.sparse.csr_array(([10.], ([0], [0])), shape=(1, 1), dtype=np.complex64)
    C = scipy.sparse.csr_array(([1.], ([0], [0])), shape=(1, 1), dtype=np.complex64)
    mats = [PETSc.Mat().createAIJ(size=K.shape, csr=(K.indptr, K.indices, K.data)) for K in (A,B,C)]

    pep.setOperators(mats)
    print('set')
    pep.solve()
    print(nconv := pep.getConverged())

    xr, xi = mats[0].createVecs()

    for i in range(nconv):
        k = pep.getEigenpair(i, xr, xi)
        error = pep.computeError(i)
        
        print("%9f%+9f j    %12g" % (k.real, k.imag, error))
        print(np.array(xr))
        
