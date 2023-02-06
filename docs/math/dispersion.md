# Dispersion

For a optical monochromatic plane wave propagating in $x_3$-direction, the $E$-field can be defined as {cite}`Liu2009`

$$
    E(x_1,x_2,x_3) = \mathcal{E}(x_1,x_2) \mathrm{e}^{\mathrm{i}(kx_3-\omega t)}
$$

## Phase velocity

The points of equal phase are where the exponent is constant, i.e.

$$
    kx_3 = \omega t + \mathrm{const}.
$$

From this we can calulate the velocity of points of constant phase as

$$
    v_p = \frac{\mathrm{d}x_3}{\mathrm{d}t} = \frac{\omega}{k}.
$$

As $k$ depends on the refractive index of the material, which is a function of the frequency, v_p depends on the frequency.
For $\frac{\mathrm{d}n}{\mathrm{d}\omega} \neq 0$, the material shows dispersion and thus leads to a phase velocity dispersion.
The case $\frac{\mathrm{d}n}{\mathrm{d}\omega} > 0$ ($\frac{\mathrm{d}n}{\mathrm{d}\lambda} < 0$) is called normal dispersion and the case $\frac{\mathrm{d}n}{\mathrm{d}\omega} < 0$ ($\frac{\mathrm{d}n}{\mathrm{d}\lambda} > 0$) is called anormalous dispersion.

## Group velocity

As usually propagating waves don't consist of a single frequency, but rather several frequencies around the carrier frequency $\omega$.
For simplifity, we consider here the case of two plane waves propagating in $x_3$-direction with equal amplitude.
The frequencies and wave propagation constants are defined as:

$$
    \omega_1 = \omega_0 + \mathrm{d}\omega \quad k_1 = k_0 + \mathrm{d}k
    \\
    \omega_2 = \omega_0 - \mathrm{d}\omega \quad k_2 = k_0 - \mathrm{d}k
$$

This way the field of the wave packet can be written as

$$
    &E
    \\ \propto
    &\mathrm{e}^{\mathrm{i}(k_1 x_3-\omega_1 t)} + \mathrm{e}^{\mathrm{i}(k_2 x_3-\omega_2 t)}
    \\ \propto
    &\cos\left(\mathrm{d}k x_3-\mathrm{d}\omega t\right) \mathrm{e}^{\mathrm{i}(k_0 x_3-\omega_0 t)}
$$

Thus, the wave has a *carrier* with the frequency $\omega_0$ and the propagation constant $k_0$ and an *envelope* which propagates as $\cos\left(\mathrm{d}k x_3-\mathrm{d}\omega t\right)$.
Similar to the phase velocity, we find for a fixed point on the envelope

$$
    \mathrm{d}kx_3 = \mathrm{d}\omega t + \mathrm{const}.
$$

This leads to the group velocity

$$
    v_g = \frac{\mathrm{d}x_3}{\mathrm{d}t} = \frac{\mathrm{d}\omega}{\mathrm{d}k}.
$$

## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```