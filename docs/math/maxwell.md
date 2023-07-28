# Background on simulations for the electromagnetic field

### Some useful relations

We here put some useful equations connecting the (angular) frequency $\omega$, the wave length $\lambda$ and the wave number $k$ in various ways (so you don't have to think about them again). The equations include also the (regular) frequency $f$ and the speed of light $c$.

$$\begin{aligned}
    \omega=2\pi f = c k  = \frac{2\pi c}{\lambda}\\
    \lambda = \frac{c}{f} = \frac{2\pi c}{\omega} = \frac{2\pi}{k} \\
    k =\frac{2\pi}{\lambda}=\frac{2\pi f}{c} = \frac{\omega}{c}
\end{aligned}$$

## Maxwell's equations 

Starting point to describe a light field are always Maxwell's equations. We here will look at Maxwell's equations in materials, which can be described by the dielectric constant. Derivations of these equations can be found in most standard textbooks on electrodynamics. 

A general way to write the Maxwell's equations is 

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{D} &= \varrho_{\text{ext}}\\
    \mathbf{\nabla}\cdot\mathbf{B} &= 0 \\
    \mathbf{\nabla}\times\mathbf{E} &= - \frac{\partial \mathbf{B}} {\partial t} \\ 
    \mathbf{\nabla}\times\mathbf{H} &= \mathbf{J}_{ext} + \frac{\partial \mathbf{D}} {\partial t} 
\end{aligned}$$ 

with the four fields

-   $\mathbf{D}$: electric displacement

-   $\mathbf{E}$: electric field

-   $\mathbf{H}$: magnetic displacement

-   $\mathbf{B}$: magnetic field

In the equations $\varrho_{\text{ext}}$ and $\mathbf{J}_{\text{ext}}$ are the external charge and current density. These are macroscopic equations, i.e., local average over microscopic quantities.

The coupling to matter can be described in the macroscopic equations via the polarisation $\mathbf{P}$ and the magnetisation $\mathbf{M}$

$$\begin{aligned}
    \mathbf{D} &=& \varepsilon_0 \mathbf{E} + \mathbf{P}\\
    \mathbf{H}&=& \frac{1}{\mu_0} \mathbf{B}- \mathbf{M}    
\end{aligned}$$ 

with $\varepsilon_0$ the electric permittivity and $\mu_0$ the magnetic permeability. The polarisation $\mathbf{P}$ describes the dipole moment per unit cell of the material

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{P} = - \varrho_{\text{int}} \notag
\end{aligned}$$ 

Now putting these back in the Maxwell equations we obtain 

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{D}= \mathbf{\nabla}\cdot\left(\varepsilon_0 \mathbf{E} + \mathbf{P} \right) 
        \qquad \Rightarrow \qquad  \mathbf{\nabla}\cdot\mathbf{E} = \frac{1}{\varepsilon_0} \varrho\, , 
\end{aligned}$$ 

which constitutes a link between the electric field $\mathbf{E}$ and all polarization effects. Here $\rho=\rho_{\text{ext}}+\rho_{\text{int}}$ is the sum of the external and internal charges. 

For most effects in light-matter interaction, it is sufficient to consider linear and isotropic materials. Then, a linear connection between polarisation and electric field can be given by the dielectric
susceptibility $\chi$ constant (for anisotropic materials $\underline{\underline{\chi}}$ is a tensor) via 

$$\mathbf{P} = \varepsilon_0 \chi \mathbf{E} \,.$$ 

With this we can write 

$$\begin{aligned}
    \mathbf{D}= \varepsilon_0 \mathbf{E} + \mathbf{P} = ( \varepsilon_0 + \varepsilon_0 \chi) \mathbf{E}  = \varepsilon_0 (1+\chi)\mathbf{E}=  \varepsilon_0 \varepsilon_r \mathbf{E}
\end{aligned}$$ 

using the dielectric displacement $\varepsilon=1+ \chi$, also called dielectric constant. The dielectric constant is often given to quantify the response of a material to an external field, hence, it
is an important quantity. Note that for the magnetic field it holds analogously

$$\mathbf{B} = \mu_0 \mu_r \mathbf{H}$$ 

with the permeablitiy $\mu$ and for non-magnetic materials (as considered here) we have $\mu=1$. We summarize 

$$\begin{aligned}
    \varepsilon = \varepsilon_0 \varepsilon_r\\
    \mu = \mu_0 \mu_r
\end{aligned}$$ 

We remind that it holds that 

$$\begin{aligned}
    c^2=\frac{1}{\mu_0\varepsilon_0}
\end{aligned}$$ 

Now we insert the description of the materials into the Maxwell equation, only two fields remain 

$$\begin{aligned}
    \varepsilon \mathbf{\nabla}\cdot\mathbf{E} &=& \varrho\\
    \mu \mathbf{\nabla}\cdot\mathbf{H} &=& 0 \\
    \mathbf{\nabla}\times\mathbf{E} &=& - \mu \frac{\partial \mathbf{H}} {\partial t} \\ 
    \mathbf{\nabla}\times\mathbf{H} &=& \mathbf{J}_{ext} + \varepsilon \frac{\partial \mathbf{E}} {\partial t} \,.
\end{aligned}$$










## Eigenvectors propagating in $x_3$-direction

Assuming no sources and currents present, Maxwell's simplifies to

$$    \begin{aligned} &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = 0 \\
    & \nabla\cdot \left(\mu\vec{\mathcal{H}}\right) = 0  \\
    & \nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t} \\
    & \nabla\times\vec{\mathcal{H}} =  \varepsilon \frac{\partial \vec{\mathcal{E}}}{\partial t}  \end{aligned} 
$$

Label:(maxwell_no_sources)

By combining the latter two equations of {eq}`maxwell_no_sources`
we get for the $\mathcal{E}$

$$
    &
    \nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right)
    =
    0

    &
    \nabla\times \left( \frac{1}{\mu}\nabla\times\vec{\mathcal{E}} \right)
    = - \varepsilon \frac{\partial^2 \vec{\mathcal{E}}}{\partial t^2}
$$ (maxwell_telegraph)

If we restrict the problem to a 2D-plane $\Omega \in \mathbb{R}^2$ like done in
{cite}`Vardapetyan2003,Vardapetyan2002,Vardapetyan2002_2`,
i.e. a plane with $\vec{x}=(x_1,x_2)$ and
assuming propagation only in $x_3$-direction with a propagation constant $\beta$,
the equations simplify for the harmonic case with a frequency of $\omega$ to:

$$
    \mathcal{E}(\vec{x},x_3,t)
    =
    (\vec{E}(\vec{x}),E_3(\vec{x}))\mathrm{e}^{i(\beta x_3 - \omega t)}

    \mathcal{H}(\vec{x},x_3,t)
    =
    (\vec{H}(\vec{x}),H_3(\vec{x}))\mathrm{e}^{i(\beta x_3 - \omega t)}
$$

Using these, the curl can be written as

$$
    \nabla \times
    =
    \begin{pmatrix}
    0 & -i \beta & \partial_y \\
    i \beta & 0 & -\partial_x \\
    -\partial_y & \partial_x & 0
    \end{pmatrix}
$$

and the derivative with respect to time becomes

$$
    \frac{\partial}{\partial t}
    = - i \omega
$$

This leads to the equations

$$
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \times \vec{E}\right)
    - \omega^2 \epsilon \vec{E}
    + \frac{\beta^2}{\mu}\vec{E}
    + \frac{i \beta}{\mu} \nabla E_3
    = 0

    &
    \nabla \cdot \left(\frac{1}{\mu} \nabla E_3\right)
    + \omega^2 \epsilon E_3
    - i \beta \nabla \cdot \left( \frac{1}{\mu} \vec{E} \right)
    = 0

    &
    \nabla \cdot \left( \epsilon \vec{E} \right)
    + i \beta \epsilon E_3
    = 0
$$

and the boundary conditions at $\partial\Omega$,
where $\vec{n}$ is the unit vector orthogonal to the boundary:

$$
    &\vec{E} \times \vec{n} = 0

    &E_3 = 0
$$

Defining

$$
    E_3^{\text{new}} = i \beta E_3
$$

converts the problem to a eigenvalue problem with the eigenvalue $\beta^2$

$$
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \times \vec{E}\right)
    - \omega^2 \epsilon \vec{E}
    + \frac{\beta^2}{\mu}\vec{E}
    + \frac{1}{\mu} \nabla E_3^{\text{new}}
    = 0

    &
    \nabla \cdot \left(\frac{1}{\mu} \nabla E_3^{\text{new}}\right)
    + \omega^2 \epsilon E_3^{\text{new}}
    + \beta^2 \nabla \cdot \left( \frac{1}{\mu} \vec{E} \right)
    = 0

    &
    \nabla \cdot \left( \epsilon \vec{E} \right)
    + \epsilon E_3^{\text{new}}
    = 0
$$

Variational problem:

$$
    &
    \left( \frac{1}{\mu} \nabla \times \vec{E}, \nabla \times \vec{F} \right)
    - \omega^2 \left( \epsilon \vec{E}, \vec{F} \right)
    + \left( \frac{1}{\mu} \nabla E_3^{\text{new}}, \vec{F} \right)
    =
    - \beta^2 \left( \frac{1}{\mu} \vec{E}, \vec{F} \right)

    &
    \left( \epsilon \vec{E}, \nabla q \right)
    -
    \left( \epsilon E_3^{\text{new}}, q \right)
    =
    0
$$

## PML

<http://www.hade.ch/docs/report_FDFD.pdf>

## Bent Waveguides

The mode profiles of bent waveguides can be calculated using the previously
derived math with an transformed effective refractive index defined as
{cite}`AzizurRahman2013`

$$
    n_{eq}(x,y)
    =
    n(x,y) \left( 1+\frac{x}{R} \right)
$$

where $R$ is the radius of curvature in $x$-direction.

See discussion on choice of R in {cite}`Masi:10`

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

## Calculating static potentials

As in the static case

$$
    \nabla\times\vec{\mathcal{E}}
    = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}
    = 0
$$

$\mathcal{E}$ can be written as

$$
    \vec{\mathcal{E}} = -\nabla \Phi
$$ (EdivPhi)

using {eq}`maxwell`, for $\Phi$ can be found that

$$
    -\nabla\cdot \left(\varepsilon \nabla \Phi\right) = \rho
$$

from which we can derive the weakform

$$
    \left(
        \varepsilon \nabla \Phi
        ,
        \nabla v
    \right)
    = \rho v
$$

which is used to calculate the potential for a given structure.
Using {eq}`EdivPhi` the electric field can be calculated from the potential.

## 2D Periodic

From {eq}`maxwell_telegraph` we get for the transverse electric field $\Psi$ {cite}`Notaros2015`

$$
    \left( \partial_x^2 + \partial_y^2 + k_0^2 n^2(x,y) \right) \Psi(x,y) = 0
$$

with

$$
    \Psi(x,y) = \mathrm{e}^{\mathrm{i}kx}\Phi(x,y)
$$

leads to

$$
    \left(
        \partial_x^2 + \partial_y^2 + \mathrm{i}2k\partial_x - k^2 + k_0^2 n^2(x,y)
    \right) \Phi(x,y)
    =
    0
$$

## Bibliography

<https://doi.org/10.1090/S0025-5718-02-01411-4>
<http://fotonica.intec.ugent.be/download/ocs131.pdf>
<https://doi.org/10.1088/0034-4885/62/3/001>

```{bibliography}
:style: unsrt
:filter: docname in docnames
```
