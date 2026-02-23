import numpy as np
from skfem import Functional


def zr(mfr, refractive_index, wavelength):
    return np.pi * mfr**2 * refractive_index / wavelength


def r_at(z, mfr, refractive_index, wavelength):
    if z == 0:
        return np.inf
    return z * (1 + z / zr(mfr, refractive_index, wavelength))


def mfr_at(mfr, z, refractive_index, wavelength):
    return mfr * np.sqrt(1 + z / zr(mfr, refractive_index, wavelength))


def e_field_gaussian(r, z, mfr, refractive_index, wavelength):
    mfr_at_z = mfr_at(mfr, z, refractive_index, wavelength)
    k = 2 * np.pi * refractive_index / wavelength

    if z != 0:
        raise Warning("Complex integration seems not to be working right")

    return (
        1
        / mfr_at_z
        / np.sqrt(np.pi / 2)
        * np.exp(-(r**2) / mfr_at_z**2)
        * np.exp(-1j * (k * z + k * r**2 / (2 * r_at(z, mfr, refractive_index, wavelength))))
    )


def overlap(basis, E_i, E_j):
    @Functional(dtype=np.complex64)
    def overlap_integral(w):
        return w["E_i"] * np.conj(w["E_j"])

    integral = overlap_integral.assemble(basis, E_i=E_i, E_j=E_j)
    norm_i = np.sqrt(overlap_integral.assemble(basis, E_i=E_i, E_j=E_i))
    norm_j = np.sqrt(overlap_integral.assemble(basis, E_i=E_j, E_j=E_j))

    return np.abs(integral / norm_i / norm_j)
