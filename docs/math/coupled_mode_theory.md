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
    \vec{E}_\nu(\vec{x})\mathrm{e}^{i(\beta x_3 - \omega t)}

    \mathcal{H}(\vec{x},x_3,t)
    =
    \sum_\nu A_\nu(z)
    \vec{H}_\nu(\vec{x})\mathrm{e}^{i(\beta x_3 - \omega t)}
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

## Bibliography

```{bibliography}
:style: unsrt
:filter: docname in docnames
```