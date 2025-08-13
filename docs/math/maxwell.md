# Introduction to optical waveguides
An *optical waveguide* is a structure which is able to confine and guide the electromagnetic field in the optical spectrum.

An optical waveguide usually consists (but not mandatorily) in at least two region of space:

- the *core*, with an higher refractive index
- the *cladding*, with a lower refractive index

```{figure} opt-wg-geom.svg
:name: fig-opt-wg-geom
:alt: typical geometry of an optical waveguide
:align: center

example of an optical waveguide geometry composed by a core of higher refractive index and a cladding of lower refractive index.
```

The specific geometry and refractive index distribution can widely vary, however most waveguides used in optics and photonics have a geometry like that represented in {numref}`fig-opt-wg-geom`

The most famous example of an optical waveguide is the *optical fiber*. Indeed in an optical fiber the electromagnetic field is confined within the core of the fiber, which has a higher refractive index than the cladding, and the electromagntic waves can flow along the fiber.Intuitively this behavior can be explained imagining that the light inside the core of a fiber is undergoes a total internal reflection when it hits the core-cladding interface trapping the light. However, a more accurate description of the waveguide properties require a description in terms of Maxwell equations.

In this page we will briefly present such an approach and we will see the physical and mathematical framework inside which femwell operate to compute the waveguides modes. This discussion do not pretend to be completely accurate and exhaustive but it is meant to serve as a tutorial for new user not used to waveguide analysis and as a reference on the approach employed by femwell for the most experienced users.

## Maxwell's equations for dielectric materials

The starting point to describe the electromagnetic field of an optical wavguide are Maxwell's equations which in their most general form they can be written as

```{math}
:label: eq-maxwell
\begin{aligned}
    \mathbf{\nabla}\cdot\mathbf{D}\left(\mathbf{r},t\right) &= \varrho_{\text{ext}}\left(\mathbf{r},t\right), \\
    \mathbf{\nabla}\cdot\mathbf{B}\left(\mathbf{r},t\right) &= 0, \\
    \mathbf{\nabla}\times\mathbf{E}\left(\mathbf{r},t\right) &= - \frac{\partial \mathbf{B}\left(\mathbf{r},t\right)} {\partial t}, \\ 
    \mathbf{\nabla}\times\mathbf{H}\left(\mathbf{r},t\right)&= \mathbf{J}_\text{ext}\left(\mathbf{r},t\right) + \frac{\partial \mathbf{D}\left(\mathbf{r},t\right)} {\partial t}, 
\end{aligned}
```

where the four fields appearing in the equations are respectively

- $\mathbf{E}\left(\mathbf{r},t\right)$ the electric field, 
- $\mathbf{D}\left(\mathbf{r},t\right)$ the electric flux density (also known as electric displacemenet),
- $\mathbf{H}\left(\mathbf{r},t\right)$ the magnetic field intensity, and
- $\mathbf{B}\left(\mathbf{r},t\right)$ the magnetic flux density, 

while $\varrho_{\text{ext}}$ and $\mathbf{J}_{\text{ext}}$ are the free charge and current density respectively. These are macroscopic equations, i.e., local averages over microscopic quantities.

The specific electrical and magnetic properties of a material are specified through the *constitutive relationship*

```{math}
:label: eq-const
\begin{aligned}
    \mathbf{D}\left(\mathbf{r},t\right) &= \varepsilon_0 \varepsilon_\text{r} (\mathbf{r}) \mathbf{E}\left(\mathbf{r},t\right)\\
    \mathbf{B}\left(\mathbf{r},t\right) &= {\mu_0} \mu_\text{r} (\mathbf{r}) \mathbf{H}\left(\mathbf{r},t\right)
\end{aligned}
```

where $\varepsilon_0$ and $\mu_0$ are the vacuum *electric permittivity* and the *magnetic permeability* respectively while $\varepsilon_\text{r}$ and $\mu_\text{r}$ are the *relative permittivity and permeability* which describe the property of the material.

Generally optical waveguides comprise *linear* and *non-magnetic* materials, so that we can assume $\mu_\text{r} = 1$ and that $\varepsilon_r$ is a function of the position $\mathbf{r}$ only and does not depend on $\mathbf{E}$. Moreover, in most of the case the material are isotropic, which implies that $\varepsilon_r$ is a scalar function of position, however in some relevant cases the core of the waveguide is realized with an anisotropic crystal. In those cases, $\varepsilon_r$ has to be regarded as a second-rank tensor.

The last assumption we do in the case of dielectric material is the absence of free charges and current, that is $\varrho_\text{ext} = 0$ and $\textbf{J}_\text{ext} = 0$

By substituting [](eq-const) into [](eq-maxwell) the Maxwell equations can be written in terms of $E$ and $H$ only.

```{math}
:label: eq-maxwell-dielectric
\begin{aligned}
    \mathbf{\nabla} \cdot \mathbf{D}\left(\mathbf{r},t\right) &= 0 \\
    \mathbf{\nabla} \cdot \mathbf{H} (\mathbf{r},t) &= 0, \\
    \mathbf{\nabla}\times\mathbf{E}(\mathbf{r},t) &= - \mu_0 \frac{\partial \mathbf{H}(\mathbf{r},t)} {\partial t}, \\ 
    \mathbf{\nabla}\times\mathbf{H}(\mathbf{r},t)&= \varepsilon_0 \varepsilon_\text{r} (\mathbf{r}) \frac{\partial \mathbf{E}(\mathbf{r},t)} {\partial t}, 
\end{aligned}
```

## Wave equation

As many users may be familiar from electromagnetic course, it is not always necessary to deal with the full set of Maxwell equations and both the electric and the magnetic fields as in [](eq-maxwell-dielectric). Indeed, if we consider 

```{math}
\begin{aligned}
    \mathbf{\nabla} \times \mathbf{\nabla} \times \mathbf{E} = \mathbf{\nabla}\times \left( - \mu_0 \frac{\partial \mathbf{H}}{\partial t}\right) = -\mu_0  \frac{\partial }{\partial t} \left(\mathbf{\nabla} \times \mathbf{H}\right) = 
        -\mu_0 \varepsilon_0 \varepsilon_\text{r} \frac{\partial^2 \mathbf{E}}{\partial t^2} 
\end{aligned}
```

Using the vector identity $\mathbf{\nabla}\times \mathbf{\nabla}\times \mathbf{E} =  \mathbf{\nabla} (\mathbf{\nabla}\cdot \mathbf{E}) - \mathbf{\nabla}^2 \mathbf{E}$ and exploiting the fact that  $\mathbf{\nabla} \times \mathbf{E} = 0$ we get to

```{math}
:label: eq-wave
\mathbf{\nabla}^2 \mathbf{E} - \frac{1}{c^2} \varepsilon_\text{r} \frac{\partial^2 \mathbf{E}}{\partial t^2} =0.
```

In this way we have obtained an equation in the $E$ field only and we have seen that this equation is indeed the wave equation.

Since in optics and photonics we want to study a system when excited with a monochromatic source at frequency $\omega$ it is convenient to study equation [](eq-wave) when the electric field oscillates sinusoidally. For this reason from now on we adopt a [phasor representation](https://en.wikipedia.org/wiki/Phasor)

```{math}
:label: eq-ansatz-phasor
\mathbf{E}(\mathbf{r}, t) = \Re \left[ \tilde{\mathbf{E}}(\mathbf{r}, \omega) e^{i- \omega t} \right]
```

By substituting the ansatz in Eq. [](eq-ansatz-phasor) into [](eq-wave), we obtain the *Helmoltz equation*

```{math}
:label: eq-helm
\mathbf{\nabla}^2 \tilde{\mathbf{E}} + \frac{\omega^2}{c^2} n^2(\textbf{r}, \omega) \tilde{\mathbf{E}} =0,
```

where we have introduced the *refractive index* $n$ defined as $\varepsilon_\text{r}(\mathbf r, \omega) = n^2 (\mathbf r, \omega)$. The phasor representation has the advantage that it contains derivatives only with respect to the space and not with the respect of time, while the frequency $\omega$ is just a parameter in the equation. Moreover the Helmoltz equation is *elliptic* making its solution through the finite element method easier to treat compared to the wave equation which is *iperbolic* \cite{QUARTERONI}. For these reason from now on we will always adopt the phasor representation as in [](eq-ansatz-phasor) and focus on the Helmoltz equation only. Then for keeping the notation simple we will always write $\mathbf{E}$ or $\mathbf{H}$ in place of $\tilde{\mathbf{E}}$ and $\tilde{\mathbf{H}}$.

## Boundary Conditions
% I would move this section to another page and maybe discuss it later or while explaining the finite element method 
<!-- It is instructive to consider the well-known case of an interface $I$ between two dielectric materials, which appear in many devices. We assume an interface between two materials called $1$ with dielectric constant $\varepsilon_1$ and $2$ with dielectric constant $\varepsilon_2$. The surface is defined by the normal vector of the interface $\mathbf{n}_{I}$ and there are no external surface charges or currents. For simplicity, we surpress the dependencies $\left(\mathbf{r},t\right)$ here. 

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

These equations are sometimes referred to as boundary conditions. It should be noted, that here boundary refers to the behaviour of the fields at the interface between two materials in contrast to boundary conditions at the edge of a simulation box.While the calculations can be done with keeping $\varepsilon(\mathbf{r})$, the boundary conditions can give a useful sanity check. They also indicate, that at interface a fine grid is required, while at areas of homogeneous materials larger grid can be chosen.  -->

## Exploiting waveguide symmetry

```{figure} opt-wg-symm.svg
:name: fig-opt-wg-symm
:alt: 3D view of a waveguide with definition of the coordinate axis
:align: center

An optical waveguide has a constant cross-section along the $z$ axis making it translational symmetric.
```

Before attempting to analyse and solve the Helmoltz equation [](eq-helm) in the specific case of optical waveguides, it is convenient to exploit the specific symmetry of these structures. Indeed, as you can see from {numref}`fig-opt-wg-symm` a waveguide has a constant cross-section for all its length and this gives the waveguide a translational symmetry that can make Eq. [](eq-helm) easier to solve.



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
