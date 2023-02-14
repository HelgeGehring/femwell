# Thermal solver

Heat transfer equation

$$
    C_p \rho \frac{\mathrm{d}T}{\mathrm{d}t}-\nabla(k\nabla T) = Q
$$

where
$C_p$ [J/kg · K] is the specific heat capacity of the material,
$\rho$ [kg/m3] is the material density,
$\frac{\mathrm{d}T}{\mathrm{d}t}$ [K/s]
is the temperature change rate of the material over time,
$k$ [W/m · K] is the thermal conductivity,
$T$ [K/m] is the temperature and
$Q$ [W/m3] is the heat per area from the source.

By introducing the thermal diffusivity

$$
    \alpha = \frac{k}{\rho C_p}
$$

The equation simplifies to

$$
    \frac{k}{\alpha} \frac{\mathrm{d}T}{\mathrm{d}t}-\nabla(k\nabla T) = Q
$$

## Steady state

At the steady state the temperature is constant,
i.e. $\frac{\mathrm{d}T}{\mathrm{d}t} = 0$ which leads to

$$
    -\nabla(k\nabla T) = Q
$$

for converting this to a weak form, we multiply the function
with a test function $v_T$ and integrate over the domain $\Omega$

$$
    \int_\Omega-\nabla(k\nabla T)v_T\mathrm{d}x = \int_\Omega Qv_T\mathrm{d}x
$$

using partial integration and choosing the boundaries condition so,
that the integral over the boundaries vanishes we get

$$
    &k\int_\Omega\nabla T\nabla v_T\mathrm{d}x = \int_\Omega Qv_T\mathrm{d}x

    \Rightarrow
    &k\nabla T\nabla v_T = Qv_T
$$

which is used for the calculation of the steady state.

## Transient state

For the description transient state,
the change of the time needs to be kept in the equation,
yielding with a similar derivation:

$$
    &
    \int_\Omega C_p \rho \frac{\mathrm{d}T}{\mathrm{d}t} v_T \mathrm{d}x
    +
    k\int_\Omega\nabla T\nabla v_T\mathrm{d}x
    =
    \int_\Omega Qv_T\mathrm{d}x

    \Rightarrow
    &
    \frac{k}{\alpha} \frac{\mathrm{d}T}{\mathrm{d}t} v_T
    +
    k\nabla T\nabla v_T = Qv_T

    \Rightarrow
    &
    \frac{k}{\alpha} \frac{\mathrm{d}T}{\mathrm{d}t} v_T
    =
    Qv_T - k\nabla T\nabla v_T
$$

which is used for the calculation of the transient behaviour.