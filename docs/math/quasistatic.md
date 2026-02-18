# Electromagnetic field simulations in the quasistatic case

The below follows the treatment of \cite{Haus1989}, Ch 3.

The quasistatic equations are derived by approximating one of the time derivatives in Maxwell's equations:
* $\frac{\partial \mathbf{B}\left(\mathbf{r},t\right)} {\partial t} \sim 0$: neglecting magnetic induction, leading to "electroquasistatic" (EQS)
* $\frac{\partial \mathbf{D}\left(\mathbf{r},t\right)} {\partial t} \sim 0$: neglecting electric displacement current, leading to "magnetoquasistatic" (MQS) 

These are justified if the condition

$$ \frac{\mu \epsilon L^2}{\tau^2} \ll 1 \rightarrow \frac{L}{v} << \tau$$

i.e. if an electromagnetic wave can propagate a characteristic length of the system in a time much shorter than the timescales $\tau$ under consideration. For instance, with $v \sim c \approx 3 \times 10^{8}$ m/s, on the chip-scale with $L \sim 1 \times 10^{-3}$ m, the timescale considered much be much larger than $\sim 1 \times 10^{-12}$ s, in practical terms below GHz frequency regimes.

Whether the electroquasistatic formulation or magnetoquasistatic formulation is used depends on which field dominates in the static limit. If the magnetic field vanishes, then the system is EQS (e.g. capacitor, non-looping conductor = resistor). Is the electric field vanishes, then the system is MQS (e.g. looping conductor = inductor). This is only easy to establish if the system is made up of "perfect conductors" and "perfect insulators", which at low frequencies regular metals and dielectrics are good approximations.


## Electroquasistatic case

In the first case, we get irrotationality (curl-free) of the E-field

$$
    \nabla\times\vec{\mathcal{E}}
    = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}
    \sim 0
$$

The lack of curl means that there is a scalar potential associated with $\mathcal{E}$:

$$
    \vec{\mathcal{E}} = -\nabla \Phi
$$ (EdivPhi)

This is the electric potential. When differences in the potential are referenced (to *e.g.* a "ground"), it is alternatively called voltage $V$.

### In insulating systems

This is explored more in-depth in \cite{Haus1989}, Ch. 4-6.

using {eq}`maxwell`, for $\Phi$ can be found that

$$
    -\nabla\cdot \left(\varepsilon \nabla \Phi\right) = \rho
$$

from which we can derive the weak form

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


### In conducting systems

This is explored more in-depth in \cite{Haus1989}, Ch. 7.

Taking the divergence of Ampère–Maxwell's law, and substituting Gauss' law, we get

$$
    \nabla \cdot \left( \nabla \times \vec{H} \right) = 0 =\nabla \cdot \left( \vec{J} + \frac{\partial \vec{D}}{\partial t} \right)
    \rightarrow 
    \nabla \cdot \vec{J} = - \frac{\partial \rho}{\partial t}
$$ (continuity)

This is also known as charge continuity. For DC currents (and without current sources), the assumption is that the net charge $\rho$ is unchanging, and hence $\nabla \cdot \vec{J} = 0$.

Many materials exhibit a linear dependence between their supported current density and the applied voltage. This is known as Ohm's law:

$$ \vec{J} = \sigma \vec{E} $$

Hence, in a purely conducting system, the electric field is the solution to

$$ -\nabla \cdot \left( \sigma \nabla \Phi \right) = 0 $$

This has a similar weak form

$$
    \left(
     \sigma \nabla \Phi
        ,
        \nabla v
    \right)
    = 0
$$

### In combined systems

When $\frac{\partial \rho}{\partial t}$ is nonzero, it can be included alongside $\vec{J}$ when taking the divergence. This leads to the displacement current $\frac{\partial \vec{D}}{\partial t}$ being added to the current:

$$
    -\nabla \cdot \left( \sigma \nabla \Phi + \frac{\partial \left( \epsilon_0 \epsilon \nabla \Phi \right)}{\partial t} \right) = 0
$$

or equivalently in Fourier space:

$$ 
    -\nabla \cdot \left( \left( \sigma + i\omega \epsilon_o \epsilon \right) \nabla \Phi \right) = 0
$$

This is exactly like the above with modified conduction. However, at DC ($\omega = 0$), in an inhomogeneous system, the current distribution is purely set by the relative conductivities. At low but nonzero frequencies, the permittivity acts as a frequency-dependent complex conductivity.


### Solving for H

Once the potential (and hence E-field) has been calculated, the magnetic field can be calculated by using the remaining equations:




## Magnetoquasistatic case

This is explored more in-depth in \cite{Haus1989}, Ch. 8-10, as well as \cite{Larson2013}, Ch. 13.

In the second case, we get

$$
    \nabla\times\vec{\mathcal{H}}
    = \varepsilon \frac{\partial \vec{\mathcal{E}}}{\partial t} + \vec{J}
    \sim \vec{J}
$$ (magnetic)

combined with the no-monopole law $\nabla \cdot \vec{B} = 0$, we see that the magnetic field is solenoidal (curl, without divergence).

The lack of divergence means that there is a vector potential associated with $\mathcal{B}$:

$$
    \vec{\mathcal{B}} = \nabla \times \vec{A}
$$ (HcurlA)

This is the magnetic vector potential. Previously in the electrostatic case we had a scalar potential, with arbitrarily-referenced potential $\Phi = 0$. Here we have a vector potential, with arbitrarily-referenced divergence $\nabla \cdot \vec{A} = 0$. This is the Coulomb gauge, chosen to simplify the below.

Taking the curl of {eq}`magnetic`, substituting in the vector potential, and using $\nabla \cdot \vec{A} = 0$:

$$
    \nabla \times \left( \mu^{-1} \nabla \times \vec{A} \right) = \vec{J}
$$

The first and easiest way to deal with this for piecewise constant permeability is to use the identity $\nabla \times \nabla \times \vec{V} = -\Delta \vec{V} + \nabla (\nabla \cdot \vec{V})$. In the Coulomb gauge, this yields a set of scalar equations for the components of $\vec A$. In 2D when $\vec{J} = J_z \hat{z}$, this is identical to the scalar Coulomb equation, solving for a scalar potential $\vec{A} = A_z \hat{z}$:

$$ -\nabla \cdot \left( \mu^{-1} \nabla \cdot A_z \right) = J_z $$

In 3D, even though this form is in principle valid, the components remain coupled through boundary conditions. Hence, the full curl weak form with enforcement of 0 divergence through Lagrange multipliers is solved directly:

$$\begin{aligned}
    \left( \mu^{-1} \nabla \times \vec{A}, \nabla \times \vec{v} \right) - \left( v, \nabla \Lambda \right) - \left( \vec{A}, \nabla q \right) = \left( \vec{J}, \vec{v} \right) \\
\end{aligned}$$

with $\vec{v} \in H_0(curl, \Omega)$ (curl-conforming),  $q \in H_0^1(\Omega)$ (scalar), and $\left( \vec{A}, \Lambda \right) \in H_0(curl, \Omega) \times H_0^1(\Omega)$

which is used to calculate the potential for a given structure.
Using {eq}`HcurlA` the magnetic field can be calculated from the vector potential.




## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```
