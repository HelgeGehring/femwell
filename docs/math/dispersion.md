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


## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```