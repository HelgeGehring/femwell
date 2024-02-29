"""
References: 

From Chrostowski, L., & Hochberg, M. (2015). Silicon Photonics Design: From Devices to Systems. Cambridge University Press. doi: 10.1017/CBO9781316084168
    Citing:
    (1) R. Soref and B. Bennett, "Electrooptical effects in silicon," in IEEE Journal of Quantum Electronics, vol. 23, no. 1, pp. 123-129, January 1987, doi: 10.1109/JQE.1987.1073206.
    (2) Reed, G. T., Mashanovich, G., Gardes, F. Y., & Thomson, D. J. (2010). Silicon optical modulators. Nature Photonics, 4(8), 518–526. doi: 10.1038/nphoton.2010.179
    (3) M. Nedeljkovic, R. Soref and G. Z. Mashanovich, "Free-Carrier Electrorefraction and Electroabsorption Modulation Predictions for Silicon Over the 1–14- $\mu\hbox{m}$ Infrared Wavelength Range," in IEEE Photonics Journal, vol. 3, no. 6, pp. 1171-1180, Dec. 2011, doi: 10.1109/JPHOT.2011.2171930.
"""

import pickle

import numpy as np
from scipy.constants import e, epsilon_0, k


def dn_carriers(wavelength: float, dN: float, dP: float) -> float:
    """Phenomenological wavelength-dependent index perturbation from free carriers in silicon.

    Use quadratic fits for wavelengths, or better-characterized models at 1550 nm and 1310 nm.

    Args:
        wavelength: (um).
        dN: excess electrons (/cm^3).
        dP: excess holes (/cm^3).

    Returns:
        dn: change in refractive index compared to intrinsic silicon.
    """
    if wavelength == 1.55:
        return -5.4 * 1e-22 * np.power(dN, 1.011) - 1.53 * 1e-18 * np.power(dP, 0.838)
    elif wavelength == 1.31:
        return -2.98 * 1e-22 * np.power(dN, 1.016) - 1.25 * 1e-18 * np.power(dP, 0.835)
    else:
        wavelength *= 1e-6
        return -3.64 * 1e-10 * wavelength**2 * dN - 3.51 * 1e-6 * wavelength**2 * np.power(
            dP, 0.8
        )


def dalpha_carriers(wavelength: float, dN: float, dP: float) -> float:
    """Phenomenological wavelength-dependent absorption perturbation from free carriers in silicon.

    Use quadratic fits for wavelengths, or better-characterized models at 1550 nm and 1310 nm.

    Args:
        wavelength: (um).
        dN: excess electrons (/cm^3).
        dP: excess holes (/cm^3).

    Returns:
        dalpha: change in absorption coefficient (/cm)
    """
    if wavelength == 1.55:
        return 8.88 * 1e-21 * dN**1.167 + 5.84 * 1e-20 * dP**1.109
    elif wavelength == 1.31:
        return 3.48 * 1e-22 * dN**1.229 + 1.02 * 1e-19 * dP**1.089
    else:
        wavelength *= 1e-6
        return 3.52 * 1e-6 * wavelength**2 * dN + 2.4 * 1e-6 * wavelength**2 * dP


def alpha_to_k(alpha, wavelength):
    """Converts absorption coefficient (/cm) to extinction coefficient (unitless), given wavelength (um)."""
    wavelength = wavelength * 1e-6  # convert to m
    alpha = alpha * 1e2  # convert to /m
    return alpha * wavelength / (4 * np.pi)


def k_to_alpha(k, wavelength):
    """Converts extinction coefficient (unitless) to absorption coefficient (/cm), given wavelength (um)."""
    wavelength = wavelength * 1e-6  # convert to m
    alpha = 4 * np.pi * k / wavelength
    return alpha * 1e-2  # convert to /cm


def k_to_alpha_dB(k, wavelength):
    """Converts extinction coefficient (unitless) to absorption coefficient (dB/cm), given wavelength (um)."""
    wavelength = wavelength * 1e-6  # convert to m
    alpha = 4 * np.pi * k / wavelength
    return 10 * np.log10(np.exp(1)) * alpha * 1e-2  # convert to /cm


# Physical constants (in cm)
q = e  # elementary charge (C)
eps = 11.68 * epsilon_0 / 100  # relative permittivity of silicon (CV−1m−1 * 1m / 100 cm)
kB = k  # Boltzmann constant (J/K)
ni = 1e10  # intrinsic carriers (cm-3)
T = 325  # Temperature (K)

# Units
um = 1e-6
cm = 1e-2
nm = 1e-9


# PN junction physics
def built_in_voltage(NA, ND):
    """Junction built-in voltage.

    Arguments:
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
    """
    return kB * T / q * np.log(NA * ND / (ni**2))


def depletion_width(NA, ND, V):
    """Depletion width.

    Arguments:
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
        V: voltage (V)
    """
    return np.sqrt(2 * eps * (NA + ND) * (built_in_voltage(NA, ND) - V) / (q * NA * ND))


def depletion_width_n_side(NA, ND, V):
    """Depletion on n-side.

    Arguments:
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
        V: voltage (V)
    """
    return depletion_width(NA, ND, V) / (1 + ND / NA)


def depletion_width_p_side(NA, ND, V):
    """Depletion on p-side.

    Arguments:
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
        V: voltage (V)
    """
    return depletion_width(NA, ND, V) / (1 + NA / ND)


def hole_concentration_depletion_approx(x, V, xpn, NA, ND):
    """Hole concentration (depletion approximation).

    Arguments:
        x: position (x, um)
        xpn: junction position (x, um)
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
        V: voltage (V)
    """
    # Acceptors on the p-side left (majority, fully ionized)
    p = np.where(x < xpn - depletion_width_p_side(NA, ND, V), NA, 0)
    # Donors on n-side right (simplest approximation to minority carriers)
    n = np.where(x > xpn + depletion_width_n_side(NA, ND, V), ni**2 / ND, 0)
    # 0 elsewhere (depletion)
    return p + n


def electron_concentration_depletion_approx(x, V, xpn, NA, ND):
    """Electron concentration (depletion approximation).

    Arguments:
        x: position (x, um)
        xpn: junction position (x, um)
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
        V: voltage (V)
    """
    # Acceptors on the p-side left (simplest approximation to minority carriers)
    p = np.where(x < xpn - depletion_width_p_side(NA, ND, V), ni**2 / NA, 0)
    # Donors on n-side right (majority, fully ionized)
    n = np.where(x > xpn + depletion_width_n_side(NA, ND, V), ND, 0)
    # 0 elsewhere (depletion)
    return p + n


def index_pn_junction(x, xpn, NA, ND, V, wavelength):
    """Refractive index of silicon around a pn junction.

    Arguments:
        x: position (um)
        xpn: junction position (um)
        NA: acceptor concentration (cm-3)
        ND: donor concentration (cm-3)
        V: voltage (V)
        wavelength: of light (um)
    """
    x_cm = x * um / cm
    xpn_cm = xpn * um / cm
    n = dn_carriers(
        wavelength,
        electron_concentration_depletion_approx(x_cm, V, xpn_cm, NA, ND),
        hole_concentration_depletion_approx(x_cm, V, xpn_cm, NA, ND),
    )
    k = alpha_to_k(
        dalpha_carriers(
            wavelength,
            electron_concentration_depletion_approx(x_cm, V, xpn_cm, NA, ND),
            hole_concentration_depletion_approx(x_cm, V, xpn_cm, NA, ND),
        ),
        wavelength,
    )
    return n + 1j * k


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    dNs = np.logspace(15, 20, 100)

    for wavelength in [1.31, 1.55]:
        dns_electrons = dn_carriers(wavelength, dNs, 0)
        dns_holes = dn_carriers(wavelength, 0, dNs)
        plt.loglog(dNs, -1 * dns_electrons, label=f"electrons, {wavelength} um")
        plt.loglog(dNs, -1 * dns_holes, label=f"holes, {wavelength} um")
    plt.xlabel("free carrier concentration (cm-3)")
    plt.ylabel("- real n")
    plt.legend()
    plt.show()

    for wavelength in [1.31, 1.55]:
        dns_electrons = k_to_alpha_dB(
            alpha_to_k(dalpha_carriers(wavelength, dNs, 0), wavelength), wavelength
        )
        dns_holes = k_to_alpha_dB(
            alpha_to_k(dalpha_carriers(wavelength, 0, dNs), wavelength), wavelength
        )
        plt.loglog(dNs, dns_electrons, label=f"electrons, {wavelength} um")
        plt.loglog(dNs, dns_holes, label=f"holes, {wavelength} um")
    plt.xlabel("free carrier concentration (cm-3)")
    plt.ylabel("absorption coeff (dB/cm)")
    plt.legend()
    plt.show()

    fig = plt.figure()

    xs = np.linspace(-1e-6, 1e-6, 1000) * 1e2  # cm

    cns = []
    cps = []

    xpn = 0  # 250E-9 * 1E2
    NA = 2e17
    ND = 2e17

    voltages = [0, -5, -10, -15]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

    for voltage, color in zip(voltages, colors):
        cns = []
        cps = []
        for x in xs:
            cns.append(electron_concentration_depletion_approx(x, voltage, xpn, NA, ND))
            cps.append(hole_concentration_depletion_approx(x, voltage, xpn, NA, ND))

        plt.semilogy(xs * 1e4, cns, color=color, label=f"{voltage} V")
        plt.semilogy(xs * 1e4, cps, color=color, linestyle="--")

    plt.legend(loc="best")

    plt.ylim([1e1, 1e19])
    plt.title("Dashed: holes, solid: electrons")
    plt.ylabel("Conc (cm$^{-3}$)")
    plt.xlabel("x (um)")

    plt.show()
