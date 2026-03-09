# Waveguide port boundary conditions

Here we introduce the math necessary for waveguide port boundary conditions, which launch a certain mode into the waveguide and absorb the reflections{cite}`Jin2015`.

## Eigenmodes
Let's start with the simplest case: a parallel-plate waveguide with a width of $b$. In this case, the fields of the modes can simply be described by

$$
h_m(y) = \sqrt{\frac{v_m}{b}}\cos\frac{m\pi y}{b},
\quad
v_m = 
\begin{cases}
    1, &m=0 \\
    2, &m \neq 0
\end{cases}
$$

and the propagation constants are

$$
\gamma_m
=
\sqrt{
    \left(
        \frac{m \pi}{b}
    \right)^2
    -
    k_0^2
}
$$

leading to the propagating modes being described as

$$
H_{z,m}(x,y) = h_m(y)\mathrm{e}^{\gamma_m x}\,.
$$

As these modes are orthogonal to each other, they fullfill the orthogonallity relation

$$
\int_0^b h_m(y)h_{m'}(y)\mathrm{d}y = \delta_{m,m'}
$$

## Field at the interface

The field at the interface can be described using these modes as

$$
H_z(x,y)
=
H_z^{inc}(x,y)
+
\sum_{m=0}^\infty
a_m h_m(y) \mathrm{e}^{\gamma_m x},
$$

where the expansion coefficients $a_m$ can be determined using

$$
a_m
=
\mathrm{e}^{-\gamma_m x_1}
\int_0^b \left[ H_z(x_1,y)-H_z^{inc}(x_1,y) \right] h_m(y) \mathrm{d}y \,.
$$

This way we get for $H_z$ the expression

$$
H_z(x,y)
=
H_z^{inc}(x,y)
+
\sum_{m=0}^\infty
\mathrm{e}^{\gamma_m (x-x_1)} h_m(y)
\int_0^b \left[ H_z(x_1,y')-H_z^{inc}(x_1,y') \right] h_m(y') \mathrm{d}y'
,
$$

and for the derivative with respect to $x$

$$
\frac{
    \partial H_z
}{
    \partial x
}
=
\frac{
    \partial H_z^{inc}
}{
    \partial x
}
+
\sum_{m=0}^\infty
\gamma_m\mathrm{e}^{\gamma_m (x-x_1)} h_m(y)
\int_0^b \left[ H_z(x_1,y')-H_z^{inc}(x_1,y') \right] h_m(y') \mathrm{d}y'
,
$$

and at the interface $x=x_1$

$$
\left.
\frac{
    \partial H_z
}{
    \partial x
}
\right|_{x=x_1}
=
\left.
\frac{
    \partial H_z^{inc}
}{
    \partial x
}
\right|_{x=x_1}
+
\sum_{m=0}^\infty
\gamma_m
\int_0^b \left[ H_z(x_1,y')-H_z^{inc}(x_1,y') \right] h_m(y') \mathrm{d}y'
,
$$

## Boundary condition

Using this, we can write it in the form of a generalized boundary condition:

$$
\frac{
    \partial H_z
}{
    \partial \vec{n}
}
+
\gamma(H_z) = q \,,
$$

where $\vec{n}$ is the vector orthogonal to the interface. the boundary operator $\gamma$ is given by

$$
\gamma(H_z) = \sum_{m=0}^\infty \gamma_m h_m(y) \int_0^b H_z(x_1,y') h_m(y') \mathrm{d}y'
$$

and $q$ is defined as

$$
q
=
\left.
\frac{
    \partial H_z^{inc}
}{
    \partial \vec{n}
}
\right|_{x=x_1}
+
\sum_{m=0}^\infty \gamma_m h_m(y) \int_0^b H_z^{inc}(x_1,y') h_m(y') \mathrm{d}y'
$$

and simplifies for single-mode incidence of mode $n$ to

$$
q = 2 \gamma_n H_0 h_n(y) \mathrm{e}^{-\gamma_n x_1} \,,
$$

where $H_0$ is the magnitude of the incident field and $n$ is the number of the incident mode.

By setting $H_0=0$ and thus setting the right-hand side of the boundary condition to zero, this kind of boundary condition can be used as an absorbing boundary condition.

## Functional defining the finite-element simulation

Adding this boundary condition to the functional defining the simulation leads to

$$
\begin{aligned}
F(H_z)
=
&\frac{1}{2}
\int_\Omega
\left[
    \frac{
        \left(
            \frac{\partial H_z}{\partial x}
        \right)^2
        +
        \left(
            \frac{\partial H_z}{\partial y}
        \right)^2
    }{\epsilon_r}
        -
        k_0^2 \mu_r H_z^2
\right]
\mathrm{d}\Omega
\\
&+
\sum_\sigma
\int_\sigma
\left[
    \frac{1}{2} H_z \gamma(H_z) - q H_z
\right]
\mathrm{d}y
\end{aligned}
$$

where $\Omega$ is the simulation domain and the $\sigma$ are the boundaries of $\Omega$, where the waveguide port condition are applied to.

As $\lim_{m\to\infty} h_m=0$, the functional converges and only needs a limited amount of summands.

## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```