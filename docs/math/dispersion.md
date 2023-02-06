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
    E
    &\propto
    \mathrm{e}^{\mathrm{i}(k_1 x_3-\omega_1 t)} + \mathrm{e}^{\mathrm{i}(k_2 x_3-\omega_2 t)}
    \\
    &\propto
    \cos\left(\mathrm{d}k x_3-\mathrm{d}\omega t\right) \mathrm{e}^{\mathrm{i}(k_0 x_3-\omega_0 t)}
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

As the energy of a wave is porportional to its field amplitude squared, the energie is concentrated in areas where the envelope is large.
Thus, the energy (and therefore also information) travels with the group velocity, while the constant-phase wavefronts travel with the phase velocity.

## Group velocity dispersion

Because the propagation constant $k$ usually has a nonlinear dependency on the frequency $\omega$, i.e. the group velocity dispersion

$$
    \frac{\mathrm{d}^2k}{\mathrm{d}\omega^2} = \frac{\mathrm{d}}{\mathrm{d}\omega}v_g^{-1} \neq 0
$$

A dimensionless coefficient for the group velocity dispersion can be defiend as

$$
    D = c \omega \frac{\mathrm{d}^2k}{\mathrm{d}\omega^2} = \frac{2\pi c^2}{\lambda} \frac{\mathrm{d}^2k}{\mathrm{d}\omega^2} \,.
$$

The group velocity dispersion describes the behaviour of the pulse envelope during propagation. It can lead to effects such as a broadening of the pulse and a time delay between pulses of different wavelengths.
For a *positive group-velocity dispersion* ($D>0$) a long-wavelength/low-frequency pulse travels faster than a short-wavelength/high-frequency pulse. Likewise, for a *negative group-velocity dispersion* ($D<0$) a short-wavelength/high-frequency pulse travels faster than a long-wavelength/low-frequency pulse. 

For describing the propagation within an optical fiber usually another definition of the dispersion coefficient is used

$$
    D_\lambda
    =
    - \frac{2\pi c}{\lambda^2} \frac{\mathrm{d}^2k}{\mathrm{d}\omega^2}
    =
    - \frac{D}{c\lambda} \,,
$$

where usually $\mathrm{\frac{ps}{km \cdot nm}}$ is used as the unit. This way, it directly described the chromatic pulse transmission delay within the fiber.

## Definition in terms of refractive index

The dispersion can also be defined in terms of the refractive index. For this the propagation constant is written as

$$
    k = \frac{\omega}{c} n(\omega)
$$

This leads to the phase velocity

$$
    v_p = \frac{c}{n} \,,
$$

the group velocity

$$
    v_g
    =
    \frac{c}{n+\omega\frac{\mathrm{d}n}{\mathrm{d}\omega}}
    =
    \frac{c}{n-\lambda\frac{\mathrm{d}n}{\mathrm{d}\lambda}}
    \,,
$$

and the group dispersion coefficient

$$
    &D(\lambda) = \lambda^2 \frac{\mathrm{d}^2n}{\mathrm{d}\lambda^2}
    \\
    &D_\lambda(\lambda) = - \frac{\lambda}{c} \frac{\mathrm{d}^2n}{\mathrm{d}\lambda^2}
    \,.
$$

In the case of nanophotonic waveguides, the refractive index are replaced by the effective refractive index of the waveguide.

## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```