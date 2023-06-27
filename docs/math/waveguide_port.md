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
H_{z,m}(x,y) = h_m(y)\exp{\gamma_m x}\,.
$$

As these modes are orthogonal to each other, they fullfill the orthogonallity relation

$$
\int_0^b h_m(y)h_{m'}(y)\mathrm{d}y = \delta_{m,m'}
$$

## Interface

The field at the interface can be described using these modes as

$$
H_z(x,y)
=
H_z^{inc}(x,y)
+
\sum_{m=0}^\inf
a_m h_m(y) \mathrm{e}^{\gamma_m x},
$$

where the expansion coefficients $a_m$ can be determined using

$$
a_m
=
\mathrm{e}^{-\gamma_m x_1}
\int_0^b \left[ H_z(x_1,y)-H_z^{inc}(x_1,y) \right] h_m(y) \mathrm{d}y \,.
$$

We obtain for the partial derivative with respect to $x$

## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```