Maxwell
-------

.. math::
    :name: maxwell

    &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = \rho

    &\nabla\cdot \left(\mu\vec{\mathcal{H}}\right) = 0

    &\nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}

    &\nabla\times\vec{\mathcal{H}} = \varepsilon\frac{\partial \vec{\mathcal{E}}}{\partial t} + \vec{J}

Without sources:

.. math::
    :name: maxwell_no_sources

    &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = 0

    &\nabla\cdot \left(\mu\vec{\mathcal{H}}\right) = 0

    &\nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}

    &\nabla\times\vec{\mathcal{H}} = \varepsilon\frac{\partial \vec{\mathcal{E}}}{\partial t}

Leads to

.. math::
    :name: maxwell_telegraph

    &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = 0

    &\nabla\times \left( \frac{1}{\mu}\nabla\times\vec{\mathcal{E}} \right)
    =
    - \varepsilon \frac{\partial^2 \vec{\mathcal{E}}}{\partial t^2}

If we restrict the problem to 2D like done in `paper <http://dx.doi.org/10.1080/02726340290084012>`_,
i.e. a plane with :math:`\vec{x}=(x_1,x_2)` and
assuming propagation only in z-direction,
the equations simplify for the harmonic case to:

.. math::
    \mathcal{E}(\vec{x},x_3,t)=(\vec{E}(\vec{x}),E_3(\vec{x}))\mathrm{e}^{i(\omega t - \beta x_3)}

    \mathcal{H}(\vec{x},x_3,t)=(\vec{H}(\vec{x}),H_3(\vec{x}))\mathrm{e}^{i(\omega t - \beta x_3)}

Leads to the equations

.. math::
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \times \vec{E}\right)
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

defining

.. math::
    E_3^{\text{new}} = i \beta E_3

converts the problem to a eigenvalue problem with the eigenvalue :math:`\beta^2`

.. math::
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \times \vec{E}\right)
    - \omega^2 \epsilon \vec{E}
    + \frac{\beta^2}{\mu}\vec{E}
    - \frac{1}{\mu} \nabla E_3^{\text{new}}
    = 0

    &
    \nabla \cdot \left(\frac{1}{\mu} \nabla E_3^{\text{new}}\right)
    + \omega^2 \epsilon E_3^{\text{new}}
    + - \beta^2 \nabla \cdot \left( \frac{1}{\mu} \vec{E} \right)
    = 0

    &
    \nabla \cdot \left( \epsilon \vec{E} \right)
    - \epsilon E_3^{\text{new}}
    = 0

Variational problem:

.. math::
    &
    \left( \frac{1}{\mu} \nabla \times \vec{E}, \nabla \times \vec{F} \right)
    - \omega^2 \left( \epsilon \vec{E}, \vec{F} \right)
    - \left( \frac{1}{\mu} \nabla E_3^{\text{new}}, \vec{F} \right)
    =
    \beta^2 \left( \frac{1}{\mu} \vec{E}, \vec{F} \right)

    &
    \left( \epsilon, \nabla q \right) + \left( \epsilon E_3^{\text{new}}, q \right)
    = 0