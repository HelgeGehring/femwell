#################
Maxwell-equations
#################

Starting with the maxwell equations:

.. math::
    :name: maxwell

    &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = \rho

    &\nabla\cdot \left(\mu\vec{\mathcal{H}}\right) = 0

    &\nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}

    &\nabla\times\vec{\mathcal{H}} = \varepsilon\frac{\partial \vec{\mathcal{E}}}{\partial t} + \vec{J}

where :math:`\mu` and :math:`\epsilon` are assumed to be element wise constant.

*************************************************
Eigenvectors propagating in :math:`x_3`-direction
*************************************************

Assuming no sources and currents present, :eq:`maxwell` simplifies to

.. math::
    :name: maxwell_no_sources

    &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = 0

    &\nabla\cdot \left(\mu\vec{\mathcal{H}}\right) = 0

    &\nabla\times\vec{\mathcal{E}} = - \mu \frac{\partial \vec{\mathcal{H}}}{\partial t}

    &\nabla\times\vec{\mathcal{H}} = \varepsilon\frac{\partial \vec{\mathcal{E}}}{\partial t}

By combining the latter two equations of :eq:`maxwell_no_sources` we get for the :math:`\mathcal{E}`

.. math::
    :name: maxwell_telegraph

    &\nabla\cdot \left(\varepsilon\vec{\mathcal{E}}\right) = 0

    &\nabla\times \left( \frac{1}{\mu}\nabla\times\vec{\mathcal{E}} \right)
    =
    - \varepsilon \frac{\partial^2 \vec{\mathcal{E}}}{\partial t^2}

If we restrict the problem to a 2D-plane :math:`\Omega \in \mathbb{R}^2` like done in :cite:`Vardapetyan2002`,
i.e. a plane with :math:`\vec{x}=(x_1,x_2)` and
assuming propagation only in :math:`x_3`-direction with a propagation constant :math:`\beta`,
the equations simplify for the harmonic case with a frequency of :math:`\omega` to:

.. math::
    \mathcal{E}(\vec{x},x_3,t)=(\vec{E}(\vec{x}),E_3(\vec{x}))\mathrm{e}^{i(\beta x_3 - \omega t)}

    \mathcal{H}(\vec{x},x_3,t)=(\vec{H}(\vec{x}),H_3(\vec{x}))\mathrm{e}^{i(\beta x_3 - \omega t)}

Using these, the curl can be written as

.. math::
    \nabla \times
    =
    \begin{pmatrix}
    0 & -i \beta & \partial_y \\
    i \beta & 0 & -\partial_x \\
    -\partial_y & \partial_x & 0
    \end{pmatrix}

and the derivative with respect to time becomes

.. math::
    \frac{\partial}{\partial t}
    =
    - i \omega

This leads to the equations

.. math::
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \times \vec{E}\right)
    - \omega^2 \epsilon \vec{E}
    + \frac{\beta^2}{\mu}\vec{E}
    + \frac{i \beta}{\mu} \nabla E_3
    = 0

    &
    \nabla \cdot \left(\frac{1}{\mu} \nabla E_3\right)
    + \omega^2 \epsilon E_3
    - i \beta \nabla \cdot \left( \frac{1}{\mu} \vec{E} \right)
    = 0

    &
    \nabla \cdot \left( \epsilon \vec{E} \right)
    + i \beta \epsilon E_3
    = 0

and the boundary conditions at :math:`\partial\Omega`,
where :math:`\vec{n}` is the unit vector orthogonal to the boundary:

.. math::
    &\vec{E} \times \vec{n} = 0

    &E_3 = 0

Defining

.. math::
    E_3^{\text{new}} = i \beta E_3

converts the problem to a eigenvalue problem with the eigenvalue :math:`\beta^2`

.. math::
    &
    \nabla \times \left(\frac{1}{\mu} \nabla \times \vec{E}\right)
    - \omega^2 \epsilon \vec{E}
    + \frac{\beta^2}{\mu}\vec{E}
    + \frac{1}{\mu} \nabla E_3^{\text{new}}
    = 0

    &
    \nabla \cdot \left(\frac{1}{\mu} \nabla E_3^{\text{new}}\right)
    + \omega^2 \epsilon E_3^{\text{new}}
    + \beta^2 \nabla \cdot \left( \frac{1}{\mu} \vec{E} \right)
    = 0

    &
    \nabla \cdot \left( \epsilon \vec{E} \right)
    + \epsilon E_3^{\text{new}}
    = 0

Variational problem:

.. math::
    &
    \left( \frac{1}{\mu} \nabla \times \vec{E}, \nabla \times \vec{F} \right)
    - \omega^2 \left( \epsilon \vec{E}, \vec{F} \right)
    + \left( \frac{1}{\mu} \nabla E_3^{\text{new}}, \vec{F} \right)
    =
    - \beta^2 \left( \frac{1}{\mu} \vec{E}, \vec{F} \right)

    &
    \left( \epsilon, \nabla q \right) - \left( \epsilon E_3^{\text{new}}, q \right)
    = 0

***
PML
***
`link <http://www.hade.ch/docs/report_FDFD.pdf>`_


***************************
TE/TM Polarization Fraction
***************************

.. math::

    \mathrm{TEfrac}
    &=
    \frac{
        \int \left| E_{x_1} \right|^2 \mathrm{d}x\mathrm{d}y
    }{
        \int \left| E_{x_1} \right|^2 + \left| E_{x_2} \right|^2 \mathrm{d}x \mathrm{d}y
    }

    \mathrm{TMfrac}
    &=
    \frac{
        \int \left| E_{x_2} \right|^2 \mathrm{d}x\mathrm{d}y
    }{
        \int \left| E_{x_1} \right|^2 + \left| E_{x_2} \right|^2 \mathrm{d}x \mathrm{d}y
    }

*********************
Loss per meter [dB/m]
*********************

.. math::
    \text{Loss at }x_3\text{ [dB]}
    &=-10 \log_{10} \frac{\left|E(x_3)\right|^2}{\left|E(x_3=0)\right|^2}
    \\
    &=-20 \log_{10} \frac{\left|E(x_3)\right|}{\left|E(x_3=0)\right|}
    \\
    &=-20 \log_{10} \mathrm{e}^{\Im\beta x_3}
    \\
    &=-20 \frac{\log_{\mathrm{e}} \mathrm{e}^{\Im\beta x_3}}{\ln 10}
    \\
    &=\frac{-20}{\ln 10} \Im\beta x_3
    \\
    \\
    \text{Loss [dB/m]}
    &=
    \frac{-20}{\ln 10} \Im\beta \, 1\mathrm{m}

**************
Effective Area
**************

As defined in :cite:p:`Agrawal2019`

.. math::
    A_{\text{eff}}
    =
    \frac{
        \left( \int \left| E \right|^2 \mathrm{d}A \right)^2
    }{
        \int \left| E \right|^4 \mathrm{d}A
    }

*******************
Overlap coefficient
*******************

.. math::
    c_{\nu\mu}
    =
    \int \mathcal{E}_\nu^* \times \mathcal{H}_\mu + \mathcal{E}_\nu \times \mathcal{H}_\mu^* \mathrm{d}A
    =
    c_{\mu\nu}^*

.. bibliography::
