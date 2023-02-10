import sys
from typing import Optional

import numpy as np
import scipy.sparse as sp
from numpy import ndarray

if "pyodide" in sys.modules:
    from scipy.sparse.base import spmatrix
else:
    from scipy.sparse import spmatrix

from skfem.utils import CondensedSystem, bmat


def mpc_symmetric(
    A: spmatrix,
    b: ndarray,
    S: Optional[ndarray] = None,
    M: Optional[ndarray] = None,
    T: Optional[spmatrix] = None,
    g: Optional[ndarray] = None,
) -> CondensedSystem:
    """Apply a multipoint constraint on the linear system.

    Parameters
    ----------
    A
    b
        The linear system to constrain.
    S
    M
    T
    g
        The constraint is of the form `x[S] = T @ x[M] + g`.

    """
    if M is None:
        M = np.array([], dtype=np.int64)
    if S is None:
        S = np.array([], dtype=np.int64)

    U = np.setdiff1d(np.arange(A.shape[0], dtype=np.int64), np.concatenate((M, S)))

    if T is None:
        T = sp.eye(len(S), len(M))
    if g is None:
        g = np.zeros(len(S))

    if T.shape[0] != len(S) or T.shape[1] != len(M) or len(g) != len(S):
        raise ValueError("Inputs to mpc have incompatible shapes.")

    B = bmat(
        [
            [
                A[U][:, U],
                A[U][:, M] + A[U][:, S] @ T,
            ],
            [
                T.T @ A[S][:, U] + A[M][:, U],
                (A[M][:, M] + T.T @ A[S][:, M] + A[M][:, S] @ T + T.T @ A[S][:, S] @ T),
            ],
        ],
        "csr",
    )

    if b.ndim == 1:
        y = np.concatenate((b[U] - A[U][:, S] @ g, b[M] - A[M][:, S] @ g))
    else:
        if np.any(np.nonzero(g)):
            raise NotImplementedError("Not yet implemented for g != 0")
        y = bmat(
            [
                [
                    b[U][:, U],
                    b[U][:, M] + b[U][:, S] @ T,
                ],
                [
                    T.T @ b[S][:, U] + b[M][:, U],
                    (b[M][:, M] + T.T @ b[S][:, M] + b[M][:, S] @ T + T.T @ b[S][:, S] @ T),
                ],
            ],
            "csr",
        )

    return (
        B,
        y,
        np.zeros(b.shape[0], dtype=B.dtype),
        (
            np.concatenate((U, M, S)),
            lambda x: np.concatenate((x, T @ x[len(U) :] + g)),
        ),
    )
