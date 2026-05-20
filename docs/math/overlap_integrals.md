# Overlap integrals

Let's write the modes of the i-th individual waveguides as $\left| i \right>_w$
and
the supermodes arising when coupling closeby-waveguides as $\left| i \right>_s$.

The overlap integral, which defines the coupling between two modes is defined as

$$
\left<\nu\,\right|\left.\mu\right>
=
\int \vec{\mathcal{E}}_\nu^* \times \vec{\mathcal{H}}_\mu
     +
     \vec{\mathcal{E}}_\nu \times \vec{\mathcal{H}}_\mu^*
\mathrm{d}A
$$

whereas the modes are normalized to
$
\left<\nu\,\right|\left.\nu\right> = 1
$

The coupling of the power from one mode to another is given by the squared overlap integral as

$$
\eta
=
P_{\text{out},\mu}/P_{\text{in},\nu}
=
\left| \left<\nu\,\right|\left.\mu\right> \right|^2
$$

As soon as we go from a single waveguide to an array of waveguides,
we do a change of basis from the basis of the individual waveguides, to the supermodes.
We do this using the operator:

$$
\sum_i \left| i \right>_s \phantom{\langle}_s\left< i \right|
$$

This leads to

$$
\left| \nu \right>_w
=
\sum_i \left| i \right>_s \phantom{\langle}_s\left< i \right| \left. \nu \right>_w
=
\sum_i c_{i\nu} \left| i \right>_s
$$
with $c_{i\nu} = \phantom{\langle}_s\left< i \right| \left. \nu \right>_w$.

The minimum transmission from the mode of a single waveguide to the mode of the same single waveguide  through this waveguide array can be calculated as (as for infinite long propagation the phase of the coefficients can be arbitrary)

$$
\left( \max_i \left|c_{\nu,i}\right|^2 - (1-\max_i \left|c_{\nu,i}\right|^2) \right)^2
=
\left( 2 \max_i \left|c_{\nu,i}\right|^2 - 1 \right)^2
$$

in case $\max_i \left|c_{\nu,i}\right|^2 > \frac{1}{2}$, otherwise, the transmission can go down to $0$.

$$
1-\left( 2 \max_i \left|c_{\nu,i}\right|^2 - 1 \right)^2
$$
