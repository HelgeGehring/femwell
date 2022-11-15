Maxwell
-------

.. math::
    :name: maxwell

    &\nabla\cdot\vec{\mathcal{E}} = \frac{\rho}{\varepsilon} \label{test}

    &\nabla\cdot\vec{\mathcal{H}} = 0

    &\nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}

    &\nabla\times\vec{\mathcal{H}} = \varepsilon\frac{\partial \vec{\mathcal{E}}}{\partial t} + \vec{J}

Without sources:

.. math::
    :name: maxwell_no_sources

    &\nabla\cdot\vec{\mathcal{E}} = 0

    &\nabla\cdot\vec{\mathcal{H}} = 0

    &\nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}

    &\nabla\times\vec{\mathcal{H}} = \varepsilon\frac{\partial \vec{\mathcal{E}}}{\partial t}

Leads to

.. math::
    :name: maxwell_telegraph

    &\nabla\cdot\vec{\mathcal{E}} = 0

    &\frac{1}{\mu} \nabla\times\nabla\times\vec{\mathcal{E}}
    =
    - \varepsilon \frac{\partial^2 \vec{\mathcal{E}}}{\partial t^2}

If we restrict the problem to 2D like done in `paper <http://dx.doi.org/10.1080/02726340290084012>`_,
i.e. a plane with :math:`\vec{x}=(x_1,x_2)` and
assuming propagation only in z-direction,
the equations simplify for the harmonic case to:

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

defining

.. math::
    E_3 = E_3^{\text{new}} = i \beta E_3

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