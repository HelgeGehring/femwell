# Electromagnetic field simulations

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

## Wave equation

From Maxwell's equation equations we can derive the wave equation. We start in the simplest case when there are no external sources, i.e., we set $\rho = 0$ and $\vec{J}=0$ and assume homogeneous $\varepsilon$ and $\mu$. Now we combine two of Maxwell's equations, namely $\vec{\nabla}\times \vec{E} = - \mu \frac{\partial \vec{H}}{\partial t}$ and
$\varepsilon \vec{\nabla}\times \vec{H} = \frac{\partial \vec{E}}{\partial t}$ in the
following way 

$$\begin{aligned}
    \vec{\nabla}\times \vec{\nabla}\times \vec{E} = \vec{\nabla}\times \left( - \mu \frac{\partial H}{\partial t}\right) = 
        -\mu \varepsilon \frac{\partial^2 \vec{E}}{\partial t^2} 
\end{aligned}$$ 

For the rotation we can use the known vector identity for the nable opertaor $\vec{\nabla}\times \vec{\nabla}\times \vec{E} =  \vec{\nabla} (\vec{\nabla}\cdot \vec{E}) - \Delta \vec{E}$.
We can further make use of Maxwell's equation that $\vec{\nabla} \vec{E} = 0$ in the case without sources. This leads us to the wave equation

$$\begin{aligned}
    \Delta \vec{E} - \frac{1}{c_{n}} \frac{\partial^2 \vec{E}}{\partial t^2} =0
\end{aligned}$$

We have introduced the important relation that the
light velocity in a medium is given by $c_n=1/\sqrt{\varepsilon\mu}$. In
vacuum, we obtain the speed of light $c_n=1/\sqrt{\varepsilon_0\mu_0}$,
while in matter the velocity is reduced by the refractive index via
$c_n=c/n$.

Solutions of the wave equation can be given in the basis of
monochromatic plane waves 

$$\begin{aligned}
    \vec{E} &=& \vec{E}_0 \cos(\vec{k} \cdot \vec{r} - \omega t + \varphi) \,, \notag 
\end{aligned}$$ 

Here, $\vec{k}$ is the wave vector and indicates the propagation direction. $E_0$ is the amplitude of the wave. Because of Maxwell's equation in free space, all waves are transversal, i.e., the
amplitude vector is perpendicular to the propagation direction $\vec{E}_0 \perp \vec{k}$. Note that this can be different in matter, in particular for nanostructured systems. The frequency is denoted by $\omega$ and there can be an additional phase $\phi$.

In many situations it is useful to write the solution as complex light
field 

$$\begin{aligned}
    \vec{E} &=& \vec{\tilde{E}}_0 e^{i(\vec{k} \cdot \vec{r}- \omega t  + \varphi)} \,
\end{aligned}$$


### Interface between piecewise constant materials 
It is instructive to consider the well-known case of an interface $I$ between two dielectric materials, which appear in many devices. We assume an interface between two materials called $1$ with dielectric constant $\varepsilon_1$ and $2$ with dielectric constant $\varepsilon_2$. The surface is defined by the normal vector of the interface $\mathbf{n}_{I}$ and there are no external surface charges or currents. For simplicity, we surpress the dependencies $\left(\mathbf{r},t\right)$ here. 

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

These equations sometimes referd to as boundary conditions. It should be noted, that here boundary refers to the behaviour of the fields at the interface between two materials in contrast to boundary conditions at the edge of a simulation box.While the calculations can be done with keeping $\varepsilon(\mathbf{r})$, the boundary conditions can give a useful sanity check. They also indicate, that at interface a fine grid is required, while at areas of homogeneous materials larger grid can be chosen. 


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
{cite}`AzizurRahman2013,shyroki2006exact,Jedidi2007,Xiao2012,Dehghannasiri2017`

$$
    n_{eq}(x,y)
    =
    n(x,y) \left( 1+\frac{x}{R} \right)
$$

where $R$ is the radius of curvature in $x$-direction.

See discussion on choice of R in {cite}`Masi:10`

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
