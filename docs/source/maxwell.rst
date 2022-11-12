Maxwell
-------

Math like in this `paper <http://dx.doi.org/10.1080/02726340290084012>`_

.. math::
    \mathcal{E}(\vec{x},x_3,t)=(\vec{E}(\vec{x}),E_3(\vec{x}))\mathrm{e}^{i(\omega t \pm \beta x_3)}

    \mathcal{H}(\vec{x},x_3,t)=(\vec{H}(\vec{x}),E_3(\vec{x}))\mathrm{e}^{i(\omega t \pm \beta x_3)}

Leads to the equations

.. math::
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \vec{E}\right)
    - \omega^2 \epsilon \vec{E}
    + \frac{\beta^2}{\mu}\vec{E}
    - \frac{i \beta}{\mu} \nabla E_3
    = 0

    &
    \nabla \cdot \left(\frac{1}{\mu} \nabla E_3\right)
    + \omega^2 \epsilon E_3
    + i \beta \nabla \cdot \left( \frac{1}{\mu} \vec{E} \right)
    = 0

    &
    \nabla \cdot \left( \epsilon \vec{E} \right)
    - i \beta \epsilon E_3
    = 0