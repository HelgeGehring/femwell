import numpy as np


def solver_dense(**kwargs):
    def solver(A, B):
        import scipy.linalg

        ks, xs = scipy.linalg.eig(A.todense(), B.todense())

        if kwargs["which"] == "LM":
            idx = np.abs(ks - kwargs["sigma"]).argsort()
            ks = ks[idx]
            xs = xs[:, idx]
        return ks, xs

    return solver


def solver_eigen_scipy_operator(**kwargs):
    """Solve generalized eigenproblem using SciPy (ARPACK).

    Returns
    -------
    EigenSolver
        A solver function that can be passed to :func:`solve`.

    """
    params = {"sigma": 10, "k": 5, "which": "LR"}
    params.update(kwargs)

    def solver(K, M, **solve_time_kwargs):
        params.update(solve_time_kwargs)
        import scipy.sparse.linalg
        from scipy.sparse.linalg import LinearOperator, eigs, factorized, splu

        M_inv = splu(M).solve
        ks, xs = eigs(
            LinearOperator(K.shape, matvec=lambda v: M_inv(K @ v), dtype=np.complex64),
            **params,
        )

        if params["which"] == "LR":
            idx = np.abs(np.real(ks)).argsort()[::-1]
            ks = ks[idx]
            xs = xs[:, idx]

        return ks, xs

    return solver


def solver_eigen_scipy_invert(**kwargs):
    """Solve generalized eigenproblem using SciPy (ARPACK).

    Returns
    -------
    EigenSolver
        A solver function that can be passed to :func:`solve`.

    """
    params = {"sigma": 10, "k": 5, "which": "LR"}
    params.update(kwargs)

    def solver(K, M, **solve_time_kwargs):
        params.update(solve_time_kwargs)
        import scipy.sparse.linalg

        A = scipy.sparse.linalg.inv(M) @ K

        ks, xs = scipy.sparse.linalg.eigs(A, **params)

        if params["which"] == "LR":
            idx = np.abs(np.real(ks)).argsort()[::-1]
            ks = ks[idx]
            xs = xs[:, idx]

        return ks, xs

    return solver


def solver_eigen_slepc(**kwargs):
    params = {"sigma": None, "k": 5, "which": "TM", "st": "SINVERT"}

    params.update(kwargs)

    def solver(K, M, **solve_time_kwargs):
        from petsc4py import PETSc
        from slepc4py import SLEPc

        which = {
            "LM": SLEPc.EPS.Which.LARGEST_MAGNITUDE,
            "SM": SLEPc.EPS.Which.SMALLEST_MAGNITUDE,
            "LR": SLEPc.EPS.Which.LARGEST_REAL,
            "SR": SLEPc.EPS.Which.SMALLEST_REAL,
            "LI": SLEPc.EPS.Which.LARGEST_IMAGINARY,
            "SI": SLEPc.EPS.Which.SMALLEST_IMAGINARY,
            "TM": SLEPc.EPS.Which.TARGET_MAGNITUDE,
            "TR": SLEPc.EPS.Which.TARGET_REAL,
            "TI": SLEPc.EPS.Which.TARGET_IMAGINARY,
            "ALL": SLEPc.EPS.Which.ALL,
        }
        st = {
            "CAYLEY": SLEPc.ST.Type.CAYLEY,
            "FILTER": SLEPc.ST.Type.FILTER,
            "PRECOND": SLEPc.ST.Type.PRECOND,
            "SHELL": SLEPc.ST.Type.SHELL,
            "SHIFT": SLEPc.ST.Type.SHIFT,
            "SINVERT": SLEPc.ST.Type.SINVERT,
        }

        params.update(solve_time_kwargs)

        K_ = PETSc.Mat().createAIJ(size=K.shape, csr=(K.indptr, K.indices, K.data))
        M_ = PETSc.Mat().createAIJ(size=M.shape, csr=(M.indptr, M.indices, M.data))

        eps = SLEPc.EPS().create()
        eps.setDimensions(params["k"])
        eps.setOperators(K_, M_)
        eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
        if params["st"]:
            eps.getST().setType(st[params["st"]])
        eps.setWhichEigenpairs(which[params["which"]])
        if params["sigma"]:
            eps.setTarget(params["sigma"])
        eps.solve()

        xr, xi = K_.getVecs()
        lams, xs = [], []
        for i in range(eps.getConverged()):
            val = eps.getEigenpair(i, xr, xi)
            lams.append(val)
            xs.append(np.array(xr) + 1j * np.array(xi))

        return np.array(lams), np.array(xs).T

    return solver


if __name__ == "__main__":
    import scipy.sparse
    from petsc4py import PETSc
    from slepc4py import SLEPc

    pep = SLEPc.PEP().create()

    A = scipy.sparse.csr_array(([24.0], ([0], [0])), shape=(1, 1), dtype=np.complex64)
    B = scipy.sparse.csr_array(([10.0], ([0], [0])), shape=(1, 1), dtype=np.complex64)
    C = scipy.sparse.csr_array(([1.0], ([0], [0])), shape=(1, 1), dtype=np.complex64)
    mats = [
        PETSc.Mat().createAIJ(size=K.shape, csr=(K.indptr, K.indices, K.data)) for K in (A, B, C)
    ]

    pep.setOperators(mats)
    print("set")
    pep.solve()
    nconv = pep.getConverged()
    print(nconv)

    xr, xi = mats[0].createVecs()

    for i in range(nconv):
        k = pep.getEigenpair(i, xr, xi)
        error = pep.computeError(i)

        print("%9f%+9f j    %12g" % (k.real, k.imag, error))
        print(np.array(xr))
