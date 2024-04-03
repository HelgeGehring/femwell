# Quantities of optical modes

## TE/TM Polarization Fraction

$$
    \mathrm{TEfrac}
    &=
    \frac{
        \int \left| E_{x_1} \right|^2 \mathrm{d}x\mathrm{d}y
    }{
        \int \left| E_{x_1} \right|^2 + \left| E_{x_2} \right|^2 \mathrm{d}x \mathrm{d}y
    }

    \mathrm{TMfrac}
    &=
    \frac{
        \int \left| E_{x_2} \right|^2 \mathrm{d}x\mathrm{d}y
    }{
        \int \left| E_{x_1} \right|^2 + \left| E_{x_2} \right|^2 \mathrm{d}x \mathrm{d}y
    }
$$

## Loss per meter [dB/m]

$$
    \text{Loss at }x_3\text{ [dB]}
    &=-10 \log_{10} \frac{\left|E(x_3)\right|^2}{\left|E(x_3=0)\right|^2}
    \\
    &=-20 \log_{10} \frac{\left|E(x_3)\right|}{\left|E(x_3=0)\right|}
    \\
    &=-20 \log_{10} \mathrm{e}^{\Im\beta x_3}
    \\
    &=-20 \frac{\log_{\mathrm{e}} \mathrm{e}^{\Im\beta x_3}}{\ln 10}
    \\
    &=\frac{-20}{\ln 10} \Im\beta x_3
    \\
    \\
    \text{Loss [dB/m]}
    &=
    \frac{-20}{\ln 10} \Im\beta \, 1\mathrm{m}
$$

## Effective Area

As defined in {cite}`Agrawal2019`

$$
    A_{\text{eff}}
    =
    \frac{
        \left( \int \left| \vec{\mathcal{E}} \right|^2 \mathrm{d}A \right)^2
    }{
        \int \left| \vec{\mathcal{E}} \right|^4 \mathrm{d}A
    }
$$

## Confinement coefficient

As defined in {cite}`Robinson2008`
(and generalized for varying refractive indices in the active area)

$$
    \Gamma
    =
    \frac{
        c \epsilon_0 \int n(\vec{x}) \left| \vec{\mathcal{E}} \right|^2 \mathrm{d}A
    }{
        \left( \int \vec{\mathcal{E}}^* \times \vec{\mathcal{H}}
        +
        \vec{\mathcal{E}} \times \vec{\mathcal{H}}^*
        \mathrm{d}A \right) / 2
    }
$$

## Overlap coefficient

$$
    c_{\nu\mu}
    =
    \frac{
        \int \vec{\mathcal{E}}_\nu^* \times \vec{\mathcal{H}}_\mu
        +
        \vec{\mathcal{E}}_\nu \times \vec{\mathcal{H}}_\mu^* \mathrm{d}A
    }{
        \prod_{i=\{\mu,\nu\}}
        \sqrt{
            \int \vec{\mathcal{E}}_i^* \times \vec{\mathcal{H}}_i
            +
            \vec{\mathcal{E}}_i \times \vec{\mathcal{H}}_i^* \mathrm{d}A
        }
    }
    =
    c_{\mu\nu}^*
$$

## Characteristic impedance

<https://ieeexplore.ieee.org/document/108320>

Power and current:

$$
    P_k
    =
    \delta_{jk}
    \int
    \left(
        \vec{\mathcal{E}}_j^* \times \vec{\mathcal{H}}_k
    \right) \cdot \hat{x}_3

    I_{zik} = \oint_{C_i} \mathcal{H} \ cdot
$$

Characteristic impedance:

$$
    P = I^T Z_c I

    Z_c = [I^{-1}]^T P I^{-1}
$$