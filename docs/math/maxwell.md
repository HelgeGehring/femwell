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

In the equations $\varrho_{\text{ext}}$ and $\mathbf{J}_{\text{ext}}$ are the external charge and current density. These are macroscopic equations, i.e., local averages over microscopic quantities.

The coupling to matter can be described in the macroscopic equations via the polarisation $\mathbf{P}$ and the magnetisation $\mathbf{M}$

$$\begin{aligned}
    \mathbf{D}\left(\mathbf{r},t\right) &=& \varepsilon_0 \mathbf{E}\left(\mathbf{r},t\right) + \mathbf{P}\left(\mathbf{r},t\right)\\
    \mathbf{H}\left(\mathbf{r},t\right) &=& \frac{1}{\mu_0} \mathbf{B}\left(\mathbf{r},t\right)- \mathbf{M}\left(\mathbf{r},t\right)    
\end{aligned}$$ 

with $\varepsilon_0$ the electric permittivity and $\mu_0$ the magnetic permeability. The polarisation $\mathbf{P}$ describes the dipole moment of the material connected to the internal charge density $\varrho_{\text{int}}(\mathbf{r})$

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{P}\left(\mathbf{r},t\right) = - \varrho_{\text{int}}\left(\mathbf{r},t\right)
\end{aligned}$$ 

Now putting these back in the Maxwell equations we obtain 

$$\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{D}\left(\mathbf{r},t\right)=& \mathbf{\nabla}\cdot\left(\varepsilon_0 \mathbf{E}\left(\mathbf{r},t\right) + \mathbf{P}\left(\mathbf{r},t\right) \right) \\
        \Rightarrow 
        \mathbf{\nabla}\cdot\mathbf{E}\left(\mathbf{r},t\right) =& \frac{1}{\varepsilon_0} \left(  \mathbf{\nabla}\cdot\mathbf{D}\left(\mathbf{r},t\right) -  \mathbf{\nabla}\cdot\mathbf{P}\left(\mathbf{r},t\right) \right)\\
        =& \frac{1}{\varepsilon_0} \left( \varrho_{\text{ext}}\left(\mathbf{r},t\right)+\varrho_{\text{int}}\left(\mathbf{r},t\right) \right) 
        = \frac{1}{\varepsilon_0} \varrho\left(\mathbf{r},t\right) 
\end{aligned}$$ 

which constitutes a link between the electric field $\mathbf{E}$ and all polarization effects. Here $\varrho=\varrho_{\text{ext}}+\varrho_{\text{int}}$ is the sum of the external and internal charges. 

For most effects in light-matter interaction, it is sufficient to consider linear materials, i.e., the relation between the polarisation and the electric field is given by a constant. For isotropic materials these constant is called *dielectric susceptibility* $\chi$ constant, while for anisotropic materials $\underline{\underline{\chi}}$ is a tensor. For the purpose of modelling dielectric materials, we assume that the polarization has the same time dependence that the electric field. We will keep the dependence on the space coordinate to account for composities of materials to bet modelling in a device, i.e., we account for spatially inhomogeneous (or structured) materials. 

$$\mathbf{P}\left(\mathbf{r},t\right) = \varepsilon_0 \chi\left(\mathbf{r}\right) \mathbf{E}\left(\mathbf{r},t\right) \,.$$ 

With this we can write 

$$\begin{aligned}
    \mathbf{D}\left(\mathbf{r},t\right)= & \varepsilon_0 \mathbf{E}\left(\mathbf{r},t\right) + \mathbf{P}\left(\mathbf{r},t\right) \\
     =& ( \varepsilon_0 + \varepsilon_0 \chi\left(\mathbf{r}\right)) \mathbf{E}\left(\mathbf{r},t\right)  \\
     =& \varepsilon_0 (1+\chi\left(\mathbf{r}\right))\mathbf{E}\left(\mathbf{r},t\right) \\
     =&  \varepsilon_0 \varepsilon_r\left(\mathbf{r}\right) \mathbf{E}\left(\mathbf{r},t\right)
\end{aligned}$$ 

using the dielectric displacement $\varepsilon=1+ \chi$, also called dielectric constant. The dielectric constant is often given to quantify the response of a material to an external field, hence, it is an important quantity. Note that for the magnetic field, assuming again isotropic, linear materials, it holds analogously

$$\mathbf{B}\left(\mathbf{r},t\right) = \mu_0 \mu_r\left(\mathbf{r}\right) \mathbf{H}\left(\mathbf{r},t\right)$$ 

with the permeablitiy $\mu_r$. Most materials are non-magnetic, such that $\mu_r=1$. We summarize 

$$\begin{aligned}
    \varepsilon\left(\mathbf{r}\right) = \varepsilon_0 \varepsilon_r\left(\mathbf{r}\right)\\
    \mu\left(\mathbf{r}\right) = \mu_0 \mu_r\left(\mathbf{r}\right)
\end{aligned}$$ 

We remind of the relation

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

From Maxwell's equation equations we can derive the wave equation. We start in the simplest case when there are no external sources, i.e., we set $\rho = 0$ and $\mathbf{J}=0$ and assume homogeneous $\varepsilon$ and $\mu$. Now we combine two of Maxwell's equations, namely $\mathbf{\nabla}\times \mathbf{E} = - \mu \frac{\partial \mathbf{H}}{\partial t}$ and
$\varepsilon \mathbf{\nabla}\times \mathbf{H} = \frac{\partial \mathbf{E}}{\partial t}$ in the
following way 

$$\begin{aligned}
    \mathbf{\nabla}\times \mathbf{\nabla}\times \mathbf{E} = \mathbf{\nabla}\times \left( - \mu \frac{\partial \mathbf{H}}{\partial t}\right) = 
        -\mu \varepsilon \frac{\partial^2 \mathbf{E}}{\partial t^2} 
\end{aligned}$$ 

For the rotation we can use the known vector identity for the nabla operator $\mathbf{\nabla}\times \mathbf{\nabla}\times \mathbf{E} =  \mathbf{\nabla} (\mathbf{\nabla}\cdot \mathbf{E}) - \Delta \mathbf{E}$.
We can further make use of Maxwell's equation that $\mathbf{\nabla}\cdot \mathbf{E} = 0$ in the case without sources. This leads us to the wave equation

$$\begin{aligned}
    \Delta \mathbf{E} - \frac{1}{c_{n}} \frac{\partial^2 \mathbf{E}}{\partial t^2} =0
\end{aligned}$$

We have introduced the important relation that the
light velocity in a medium is given by $c_n=1/\sqrt{\varepsilon\mu}$. In
vacuum, we obtain the speed of light $c=1/\sqrt{\varepsilon_0\mu_0}$,
while in matter the velocity is reduced by the refractive index via
$c_n=c/n$ with the refractive index $n=\sqrt{\varepsilon_r \mu_r}$.

Solutions of the wave equation can be given in the basis of
monochromatic plane waves 

$$\begin{aligned}
    \mathbf{E} &=& \mathbf{E}_0 \cos(\mathbf{k} \cdot \mathbf{r} - \omega t + \varphi) \,. 
\end{aligned}$$ 

Here, $\mathbf{k}$ is the wave vector and indicates the propagation direction. $E_0$ is the amplitude of the wave. Because of Maxwell's equation in free space, all waves are transversal, i.e., the
amplitude vector is perpendicular to the propagation direction $\mathbf{E}_0 \perp \mathbf{k}$. Note that this can be different in matter, in particular for nanostructured systems. The frequency is denoted by $\omega$ and there can be an additional phase $\varphi$.

In many situations it is useful to write the solution as complex light
field 

$$\begin{aligned}
    \mathbf{E} &=& \mathbf{\tilde{E}}_0 e^{i(\mathbf{k} \cdot \mathbf{r}- \omega t  + \varphi)} \,
\end{aligned}$$


## Interface between piecewise constant materials 
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

These equations are sometimes referred to as boundary conditions. It should be noted, that here boundary refers to the behaviour of the fields at the interface between two materials in contrast to boundary conditions at the edge of a simulation box.While the calculations can be done with keeping $\varepsilon(\mathbf{r})$, the boundary conditions can give a useful sanity check. They also indicate, that at interface a fine grid is required, while at areas of homogeneous materials larger grid can be chosen. 


## Reduction to two dimensions

Often structures are constructed in such a way, that they are strongly patterned in a plane (e.g. the $xy$-plane) and are uniform in the direction perpendicular to the plane. This already holds for the simple
example of a waveguide, which has interface within the plane, but is (almost) infinitely extended in the perpendicular direction.

For simulations, just solving the equation within the plane can reduce the problem greatly, also from the computational point of view. The equation also reduce to a much simpler version, which we derive here.

We will discuss only the electric field here. Equations for the magnetic field can be derived in a similar way. We will make a plane wave ansatz assuming a monochromatic wave with frequency $\omega$ via

$$\mathbf{E}(\mathbf{r},t) = \mathbf{E}(\mathbf{r}) e^{i(\mathbf{k} \cdot \mathbf{r}-\omega t)}$$

Now we assume that the wave propagates along the $z$-direction and, thus, $\mathbf{k} \to  k_z $. The electric field is only structured in the $x,y$-plane and depends only on $\mathbf{r}_{\perp }=(x,y) $ while in the $z$- direction the field is homogeneous in this ansatz. We separate the amplitude into the amplitude within the $xy$-plane $\mathbf{E}_{ \perp }=(E_x,E_y,0)$ and along the propagation direction $E_z$, that is we take $    \mathbf{E}(\mathbf{r})  \to  \mathbf{E} ( \mathbf{r}_{ \perp }) = \left(\mathbf{E}_{ \perp }(\mathbf{r}_{ \perp } ),E_z( \mathbf{r}_{\perp})\right) $. Note that $\mathbf{E}$ is a 2D vector in the xy plane, while $E_z$ is a scalar and that both amplitudes depend on $\mathbf{r}_{ \perp }$. 

$$
\mathbf{E}(\mathbf{r},t) \to \left(\mathbf{E}_{ \perp }(\mathbf{r}_{ \perp }),E_z(\mathbf{r}_{ \perp })\right) e^{i(k_z z- \omega t )}
$$

With these assumption we can reduce the Maxwell and wave equation respectively. For the derivations in these equations appearing we have

$$ \begin{aligned}
\frac{\partial}{\partial z}\mathbf{E}(\mathbf{r},t) = i k_z \mathbf{E}(\mathbf{r},t) \\
	\frac{\partial}{\partial t} \mathbf{E}(\mathbf{r},t) = -i \omega \mathbf{E}(\mathbf{r},t)
\end{aligned}$$

We enter this into the the wave equation
$$
\mathbf{\nabla}\times  \left( \frac{1}{\mu(\mathbf{r})} \mathbf{\nabla} \times \mathbf{E}(\mathbf{r},t)\right)  - \frac{\partial^2}{\partial t^2} \mathbf{E}(\mathbf{r},t) = 0 \,.
$$
Removing the time dependence of the fields, only the field
amplitudes enter and we obtain two equations 

$$\begin{aligned}
     \Rightarrow & \mathbf{\nabla} \times \frac{1}{\mu(\mathbf{r}_{ \perp } ) } \mathbf{\nabla} \times \mathbf{E}_{ \perp }(\mathbf{r}_{ \perp } )  
		- \left[ \omega^2 \varepsilon(\mathbf{r}_{ \perp }) +\frac{k_z^2}{\mu(\mathbf{r}_{ \perp })} \right] 
          \mathbf{E}_{\perp}(\mathbf{r}_{ \perp }) - i \frac{k_z}{\mu(\mathbf{r}_{ \perp })} \nabla E_z(\mathbf{r}_{ \perp })= 0 \\ 
     &\nabla \cdot \left( \varepsilon(\mathbf{r}_{\perp}) \nabla \cdot E_z (\mathbf{r}_{\perp}) \right) + \omega^2 \varepsilon(\mathbf{r}_{\perp})  E_z (\mathbf{r}_{\perp}) - i k_z \nabla \cdot \left( \frac{1}{\mu(\mathbf{r}_{\perp})} \mathbf{E}_{\perp}(\mathbf{r}_{\perp}) \right) = 0 
 \end{aligned} 
$$

Note that the top equation is a 2D equation within the $(x,y)$ plane, while the bottom equation is just a scalar equation.

We can proceed analogously for the Maxwell's equation 

$$
    \mathbf{\nabla} \left( \varepsilon(\mathbf{r})\mathbf{E}(\mathbf{r},t)  \right) =0
$$

and remove the time dependence from the equations simplifying it to

$$
     \Rightarrow  \mathbf{\nabla} \left( \varepsilon(\mathbf{r}_{ \perp }) \mathbf{E}_{\perp}(\mathbf{r}_{ \perp })\right)
		+ i k_z \varepsilon(\mathbf{r}_{ \perp }) E_z(\mathbf{r}_{ \perp })=0 
$$

Again this is a scalar equation. Now we have four equations for three free parameters $E_x,E_y,E_z$, hence our system is overestimated. Hence, we can skip one equation, which we chose to be the scalar equation from the wave equation. 
Defining

$$
    E_3 = i k_z E_z
$$
and dropping the indices $\perp,z$ we obtain the following two equations

$$\begin{aligned}
    &
    \mathbf{\nabla} \times \left(\frac{1}{\mu(\mathbf{r})} \mathbf{\nabla} \times \mathbf{E}(\mathbf{r})\right)
    - \omega^2 \varepsilon(\mathbf{r}) \mathbf{E}(\mathbf{r})
    - \frac{1}{\mu(\mathbf{r})} \mathbf{\nabla} E_3(\mathbf{r})
    =- \frac{k^2}{\mu(\mathbf{r})}\mathbf{E}(\mathbf{r})\\
    &
    \mathbf{\nabla} \cdot \left( \epsilon(\mathbf{r}) \mathbf{E}(\mathbf{r}) \right)
    + \epsilon E_3(\mathbf{r})
    = 0 \end{aligned}
$$
Now we have reduced the system to two dimensions by excuting the derivatives. 

### Variational eigenvalue problem
In FEM simulations, we solve these equations by a variational ansatz. That means, we take test functions, that are defined on the same space that the wave functions, and search for the optimum. To efficiently solve this, we rewrite the equations from above into an eigenvalue problem. 

For this, we define the test functions $\mathbf{F}$ and $q$ and the inner product 
$$
 \left( \mathbf{F}, \mathbf{G}\right) = \int_\Omega dr^2 \mathbf{F}\cdot \mathbf{G}^*
$$
with $g^*$ being the complex conjugated and $\Omega$ is the area. We now use this inner product to convert the equations (we drop the $\mathbf{r}$ dependence here for clarity) to

$$\begin{aligned}
    &
    \left( \frac{1}{\mu} \mathbf{\nabla} \times \mathbf{E}, \mathbf{\nabla} \times \mathbf{F} \right)
    - \omega^2 \left( \varepsilon \mathbf{E}, \mathbf{F} \right)
    + \left( \frac{1}{\mu} \mathbf{\nabla} E_3, \mathbf{F} \right)
    =
    - k^2 \left( \frac{1}{\mu} \mathbf{E}, \mathbf{F} \right)\\
    &
    \left( \epsilon \mathbf{E}, \mathbf{\nabla} q \right)
    -
    \left( \varepsilon E_3, q \right)
    =
    0  \end{aligned}
$$

This an eigenvalue-equation for $k$ and under certain assumption, it can be shown that $k=0$ is not an eigenvalues as well and the solution is stable for $\omega \to 0$. 

This is the equation we solve in the FEM.


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
