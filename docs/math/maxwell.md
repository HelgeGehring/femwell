# Background on simulations for the electromagnetic field

### Some useful relations

We here put some useful equations connecting the (angular) frequency $\omega$, the wave length $\lambda$ and the wave number $k$ in various ways (so you don't have to think about them again). The equations include also the (regular) frequency $f$ and the speed of light $c$.

$$\begin{aligned}
    \omega=2\pi f = c k  = \frac{2\pi c}{\lambda}\\
    \lambda = \frac{c}{f} = \frac{2\pi c}{\omega} = \frac{2\pi}{k} \\
    k =\frac{2\pi}{\lambda}=\frac{2\pi f}{c} = \frac{\omega}{c}
\end{aligned}$$

## Maxwell's equations for dielectric materials

Starting point to describe a light field are always Maxwell's equations. We here will look at Maxwell's equations in materials, which can be described by the dielectric constant. Derivations of these equations can be found in most standard textbooks on electrodynamics. 

A general way to write the Maxwell's equations is 

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{D}\left(\mathbf{r},t\right) &= \varrho_{\text{ext}}\left(\mathbf{r},t\right)\\
    \mathbf{\nabla}\cdot\mathbf{B}\left(\mathbf{r},t\right) &= 0 \\
    \mathbf{\nabla}\times\mathbf{E}\left(\mathbf{r},t\right) &= - \frac{\partial \mathbf{B}\left(\mathbf{r},t\right)} {\partial t} \\ 
    \mathbf{\nabla}\times\mathbf{H}\left(\mathbf{r},t\right)&= \mathbf{J}_{ext}\left(\mathbf{r},t\right) + \frac{\partial \mathbf{D}\left(\mathbf{r},t\right)} {\partial t} 
\end{aligned}$$ 

with the four fields (note that a bold letter indicates a vector)

-   $\mathbf{D}\left(\mathbf{r},t\right)$: electric displacement

-   $\mathbf{E}\left(\mathbf{r},t\right)$: electric field

-   $\mathbf{H}\left(\mathbf{r},t\right)$: magnetic displacement

-   $\mathbf{B}\left(\mathbf{r},t\right)$: magnetic field

In the equations $\varrho_{\text{ext}}$ and $\mathbf{J}_{\text{ext}}$ are the external charge and current density. These are macroscopic equations, i.e., local average over microscopic quantities.

The coupling to matter can be described in the macroscopic equations via the polarisation $\mathbf{P}$ and the magnetisation $\mathbf{M}$

$$\begin{aligned}
    \mathbf{D}\left(\mathbf{r},t\right) &=& \varepsilon_0 \mathbf{E}\left(\mathbf{r},t\right) + \mathbf{P}\left(\mathbf{r},t\right)\\
    \mathbf{H}\left(\mathbf{r},t\right) &=& \frac{1}{\mu_0} \mathbf{B}\left(\mathbf{r},t\right)- \mathbf{M}\left(\mathbf{r},t\right)    
\end{aligned}$$ 

with $\varepsilon_0$ the electric permittivity and $\mu_0$ the magnetic permeability. The polarisation $\mathbf{P}$ describes the dipole moment per unit cell of the material

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{P}\left(\mathbf{r},t\right) = - \varrho_{\text{int}}\left(\mathbf{r},t\right) \notag
\end{aligned}$$ 

Now putting these back in the Maxwell equations we obtain 

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{D}\left(\mathbf{r},t\right)=& \mathbf{\nabla}\cdot\left(\varepsilon_0 \mathbf{E}\left(\mathbf{r},t\right) + \mathbf{P}\left(\mathbf{r},t\right) \right) \\
        \Rightarrow 
        \mathbf{\nabla}\cdot\mathbf{E}\left(\mathbf{r},t\right) =& \frac{1}{\varepsilon_0} \left(  \mathbf{\nabla}\cdot\mathbf{D}\left(\mathbf{r},t\right) -  \mathbf{\nabla}\cdot\mathbf{P}\left(\mathbf{r},t\right) \right)
        = \frac{1}{\varepsilon_0} \left( \rho_{\text{ext}}\left(\mathbf{r},t\right)+\rho_{\text{int}}\left(\mathbf{r},t\right) \right) 
        = \frac{1}{\varepsilon_0} \rho\left(\mathbf{r},t\right) 
\end{aligned}$$ 

which constitutes a link between the electric field $\mathbf{E}$ and all polarization effects. Here $\rho=\rho_{\text{ext}}+\rho_{\text{int}}$ is the sum of the external and internal charges. 

For most effects in light-matter interaction, it is sufficient to consider linear materials, i.e., the relation between the polarisation and the electric field is given by a constant. For isotropic materials these constant is called *dielectric susceptibility* $\chi$ constant, while for anisotropic materials $\underline{\underline{\chi}}$ is a tensor. For the purpose of modelling dielectric materials, we assume that the polarization has the same time dependence that the electric field. We will keep the dependence on the space coordinate to account for composities of materials to bet modelling in a device, i.e., we account for inhomogeneous materials. 

$$\mathbf{P}\left(\mathbf{r},t\right) = \varepsilon_0 \chi\left(\mathbf{r}\right) \mathbf{E}\left(\mathbf{r},t\right) \,.$$ 

With this we can write 

$$\begin{aligned}
    \mathbf{D}\left(\mathbf{r},t\right)= \varepsilon_0 \mathbf{E}\left(\mathbf{r},t\right) + \mathbf{P}\left(\mathbf{r},t\right) = ( \varepsilon_0 + \varepsilon_0 \chi\left(\mathbf{r}\right)) \mathbf{E}\left(\mathbf{r},t\right)  = \varepsilon_0 (1+\chi\left(\mathbf{r}\right))\mathbf{E}\left(\mathbf{r},t\right)=  \varepsilon_0 \varepsilon_r\left(\mathbf{r}\right) \mathbf{E}\left(\mathbf{r},t\right)
\end{aligned}$$ 

using the dielectric displacement $\varepsilon=1+ \chi$, also called dielectric constant. The dielectric constant is often given to quantify the response of a material to an external field, hence, it is an important quantity. Note that for the magnetic field, assuming again isotropic, linear materials, it holds analogously

$$\mathbf{B}\left(\mathbf{r},t\right) = \mu_0 \mu_r\left(\mathbf{r}\right) \mathbf{H}\left(\mathbf{r},t\right)$$ 

with the permeablitiy $\mu_r$. For materials are non-magnetic, such that we have $\mu_R=1$. We summarize 

$$\begin{aligned}
    \varepsilon\left(\mathbf{r}\right) = \varepsilon_0 \varepsilon_r\left(\mathbf{r}\right)\\
    \mu\left(\mathbf{r}\right) = \mu_0 \mu_r\left(\mathbf{r}\right)
\end{aligned}$$ 

We remind that it holds that 

$$\begin{aligned}
    c^2=\frac{1}{\mu_0\varepsilon_0}
\end{aligned}$$ 

Now we insert the description of the materials into the Maxwell equation, only two fields remain.

$$\begin{aligned}
    \mathbf{\nabla}\cdot \left( \varepsilon\left(\mathbf{r}\right)  \mathbf{E}\left(\mathbf{r},t\right) \right) &=& \varrho \left(\mathbf{r},t\right)\\
     \mathbf{\nabla}\cdot \left( \mu \left(\mathbf{r}\right) \mathbf{H}\left(\mathbf{r},t\right) \right) &=& 0 \\
    \mathbf{\nabla}\times\mathbf{E}\left(\mathbf{r},t\right) &=& - \mu \frac{\partial \mathbf{H}\left(\mathbf{r},t\right)} {\partial t} \\ 
    \mathbf{\nabla}\times\mathbf{H}\left(\mathbf{r},t\right) &=& \mathbf{J}_{\text{ext}}\left(\mathbf{r},t\right) + \varepsilon \frac{\partial \mathbf{E}\left(\mathbf{r},t\right)} {\partial t} \,.
\end{aligned}$$

### Piecewise constant materials and interface conditions
It is instructive to condsider the well-known case of an interface $I$ between two dielectric materials, which appear in many devices. We assume an interface between two materials called $1$ with dielectric constant $\varepsilon_1$ and $2$ with dielectric constant $\varepsilon_2$. The surface is defined by the normal vector of the interface $\mathbf{n}_{I}$ and there are no external surface charges or currents. For simplicity, we surpress the dependencies $\left(\mathbf{r},t\right)$ here. 

All fields can then be split into the component parallel to the interface (hence perpendicular to the normal vector) and perpendicular to the interface (hence parallel to the normal vector). For example we consider the electric field: Define the normalized field vector $\hat{\mathbf{E}}=\mathbf{E}/E$ we split it into

$$\begin{aligned}
    \mathbf{E}=& \left(\mathbf{n}_I \times \hat{\mathbf{E}} \right) \mathbf{E} +  \left(\mathbf{n}_I \cdot \hat{\mathbf{E}} \right) \mathbf{E}\\
     =&\mathbf{E}^{\parallel}_I + \mathbf{E}^{\perp}_I
\end{aligned}$$
 
Maxwell's equation now impose continuity conditions at surfaces:

tangential component of $\mathbf{E}$ is continuous 

$$ \mathbf{E}^{\parallel}_1 = \mathbf{E}^{\parallel}_2 $$
    
normal component of $\mathbf{D}$ is continuous
 
$$\mathbf{D}^{\perp}_1 = \mathbf{D}^{\perp}_2  \quad \Leftrightarrow \quad  
                \varepsilon_1 \mathbf{E}^{\perp}_1 = \varepsilon_2 \mathbf{E}^{\perp}_2$$
   
tangential component of $\mathbf{H}$ is continuous

$$\mathbf{H}^{\parallel}_1 =  \mathbf{H}^{\parallel}_2$$
  
normal component of $\mathbf{B}$ is continuous
   
$$\mathbf{B}^{\perp}_1 =  \mathbf{B}^{\perp}_2  \quad \Leftrightarrow \quad   
                \frac{1}{\mu_1}\mathbf{H}^{\perp}_1 =  \frac{1}{\mu_2}\mathbf{H}^{\perp}_2 $$

In these equations we have used $\mathbf{D}=\varepsilon_0\varepsilon\mathbf{E}$ and $\mathbf{H}=\frac{1}{\mu_0}\mathbf{B} $.

While the calculations can be done with keeping $\varepsilon(\mathbf{r})$, the interface conditions can give a useful sanity check. They also indicate, that at interfaces a fine grid is required, while at areas of homogeneous materials larger grid can be chosen. 


## Eigenvectors propagating in $x_3$-direction

Assuming no sources and currents present, Maxwell's simplifies to

$$    \begin{aligned} &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = 0 \\
    & \nabla\cdot \left(\mu\vec{\mathcal{H}}\right) = 0  \\
    & \nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t} \\
    & \nabla\times\vec{\mathcal{H}} =  \varepsilon \frac{\partial \vec{\mathcal{E}}}{\partial t}  \end{aligned} 
$$(maxwell_no_sources)

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
