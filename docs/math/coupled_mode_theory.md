# Coupled Mode Theory

<https://www.fiberoptics4sale.com/blogs/wave-optics/two-mode-coupling>
<https://www.fiberoptics4sale.com/blogs/wave-optics/coupled-mode-theory>
<http://home.iitj.ac.in/~k.r.hiremath/research/thesis.pdf>
<http://home.iitj.ac.in/~k.r.hiremath/research/stoffer_hiremath_I_3.pdf>

## Coupling between different modes of the same waveguide

Based on {cite}`Liu2009`

$$
    \mathcal{E}(\vec{x},x_3,t)
    =
    \sum_\nu A_\nu(z)
    \vec{E}_\nu(\vec{x})\mathrm{e}^{i\beta x_3}

    \mathcal{H}(\vec{x},x_3,t)
    =
    \sum_\nu A_\nu(z)
    \vec{H}_\nu(\vec{x})\mathrm{e}^{i\beta x_3}
$$

$$  
    &
    \nabla\times\vec{\mathcal{E}}
    =
    i \mu_0 \omega \mathcal{H}

    &
    \nabla\times\vec{\mathcal{H}}
    =
    -i \omega \varepsilon \mathcal{E} - i \omega \Delta \mathcal{P}
$$

$$
    \nabla \left( \mathcal{E}_1 \times \mathcal{H}_2^* + \mathcal{H}_2^* \times \mathcal{E}_1 \right)
    =
    - i \omega \left( \mathcal{E}_1 \cdot \Delta \mathcal{P}_2^* + \mathcal{E}_2^* \cdot \Delta \mathcal{P}_1 \right)
$$

$$
    \sum_\nu \frac{\mathrm{d}}{\mathrm{d}z} A_\nu(z) \mathrm{e}^{i(\beta_\nu-\beta_\mu)z}
    \int_\Omega
        \left( E_\nu \times H_\mu^* + E_\mu^* \times H_\nu \right) \cdot \hat{x_3}
    \mathrm{d}A
    =
    - i \omega \mathrm{e}^{-i\beta_\mu z}
    \int_\Omega
        E_\mu^* \cdot \nabla P
    \mathrm{d}A
$$

$$
    \frac{\mathrm{d} A_\nu(z)}{\mathrm{d}z}
    =
    - i \omega \mathrm{e}^{-i\beta_\nu z}
    \int_\Omega
        E_\nu^* \cdot \nabla P
    \mathrm{d}A
$$

$$
    \Delta \mathcal{P}
    =
    \Delta \varepsilon \mathcal{E}
    =
    \Delta \varepsilon
    \sum_\nu A_\nu(z)
    \vec{E}_\nu(\vec{x})\mathrm{e}^{i\beta_\nu x_3}
$$

$$
    \pm
    \frac{\mathrm{d} A_\nu(z)}{\mathrm{d}z}
    =
    \sum_\nu A_\nu(z)
    i \kappa_{\nu\mu} A_\mu \mathrm{e}^{i(\beta_\mu-\beta_\nu)z}
$$

$$
    \kappa_{\nu\mu}
    =
    \omega
    \int_\Omega
        E_\nu^* \cdot \Delta \varepsilon \cdot E_\mu P
    \mathrm{d}A
$$

For lossles, i.e. $ \varepsilon = \varepsilon^* $

$$
    \kappa_{\nu\mu} = \kappa_{\mu\nu}^*
$$


## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```