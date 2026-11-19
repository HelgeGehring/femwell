"""Tests of Mode.calculate_voltage, Mode.calculate_current and Mode.calculate_impedance."""

import dataclasses
from collections import OrderedDict

import numpy as np
import scipy.constants
import shapely
import shapely.ops
from shapely.geometry import LineString, box
from skfem import (
    Basis,
    ElementTriN1,
    ElementTriN2,
    ElementTriP0,
    ElementTriP1,
    ElementTriP2,
    Functional,
    MeshTri,
)
from skfem.helpers import inner
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import Mode, compute_modes
from femwell.mesh import mesh_from_OrderedDict


def _analytic_mode(order):
    """Mode on a square mesh with polynomial fields, for which all integrals are known.

    E = (1+0.5j)*(1, 2) for order 1 and (1+0.5j)*(y+1, x+2) for order 2, i.e. E = -grad(phi) with
    phi = -(1+0.5j)*(x+2y) and phi = -(1+0.5j)*(xy+x+2y).
    H = (-y/2+0.3, x/2-0.1), which has curl_z(H) = 1.
    """
    mesh = MeshTri.init_tensor(np.linspace(-1, 1, 21), np.linspace(-1, 1, 21)).with_subdomains(
        {"cond": lambda x: (abs(x[0]) < 0.3) & (abs(x[1]) < 0.3)}
    )
    if order == 1:
        basis = Basis(mesh, ElementTriN1() * ElementTriP1())
    else:
        basis = Basis(mesh, ElementTriN2() * ElementTriP2())
    basis_t = basis.split_bases()[0]

    def field(function):
        values = basis.zeros(dtype=complex)
        values[basis.split_indices()[0]] = basis_t.project(function)
        return values

    if order == 1:
        E = (1 + 0.5j) * field(lambda x: np.array([1.0 + 0 * x[0], 2.0 + 0 * x[0]]))
    else:
        E = (1 + 0.5j) * field(lambda x: np.array([x[1] + 1.0, x[0] + 2.0]))
    H = field(lambda x: np.array([-x[1] / 2 + 0.3, x[0] / 2 - 0.1]))

    basis_epsilon_r = basis.with_element(ElementTriP0())
    mode = Mode(
        frequency=1e9,
        k=1.0,
        basis_epsilon_r=basis_epsilon_r,
        epsilon_r=basis_epsilon_r.ones(),
        basis=basis,
        E=E,
        H=H,
    )
    return mode


def _potential(order, point):
    x, y = point
    if order == 1:
        return -(1 + 0.5j) * (x + 2 * y)
    return -(1 + 0.5j) * (x * y + x + 2 * y)


def _flux(mode):
    """int (E_t x H_t*) . z dA, with complex128."""

    @Functional(dtype=complex)
    def flux(w):
        return w["E"][0][0] * np.conj(w["H"][0][1]) - w["E"][0][1] * np.conj(w["H"][0][0])

    return flux.assemble(
        mode.basis, E=mode.basis.interpolate(mode.E), H=mode.basis.interpolate(mode.H)
    )


def _raises(error, function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except error:
        return
    raise AssertionError(f"{error.__name__} not raised")


def test_voltage_exact():
    paths = (
        ((-0.83, -0.41), (0.77, 0.93)),
        ((-1, 0), (1, 0)),  # along mesh edges
        ((-0.9, -0.9), (0.9, 0.9)),
    )
    for order in (1, 2):
        mode = _analytic_mode(order)
        for start, end in paths:
            exact = _potential(order, start) - _potential(order, end)
            voltage = mode.calculate_voltage(start, end)
            np.testing.assert_allclose(voltage, exact, rtol=0, atol=1e-12)
            reverse = mode.calculate_voltage(end, start)
            np.testing.assert_allclose(reverse, -voltage, rtol=0, atol=1e-12)


def test_voltage_errors():
    mode = _analytic_mode(1)
    _raises(ValueError, mode.calculate_voltage, (0.1, 0.1), (0.1, 0.1))
    _raises(ValueError, mode.calculate_voltage, (0, 0), (1.5, 0))


def test_current_exact():
    mode = _analytic_mode(2)
    # counter-clockwise around the 0.6 x 0.6 square: area * curl_z(H) = 0.36
    np.testing.assert_allclose(mode.calculate_current(elements="cond"), 0.36, rtol=0, atol=1e-12)
    inside = lambda x: (abs(x[0]) < 0.3) & (abs(x[1]) < 0.3)  # noqa: E731
    np.testing.assert_allclose(mode.calculate_current(elements=inside), 0.36, rtol=0, atol=1e-12)
    # the +z current beyond the outer boundary is minus the current in the domain
    np.testing.assert_allclose(
        mode.calculate_current(facets=mode.basis.mesh.boundary_facets()),
        -4,
        rtol=0,
        atol=1e-12,
    )


def test_current_errors():
    mode = _analytic_mode(1)
    mesh = mode.basis.mesh
    _raises(ValueError, mode.calculate_current, elements="cond", facets=mesh.boundary_facets())
    _raises(ValueError, mode.calculate_current)
    _raises(ValueError, mode.calculate_current, facets=mesh.facets_around("cond"))


def test_impedance_formulas():
    mode = _analytic_mode(2)
    start, end = (-0.5, -0.4), (0.6, 0.7)
    voltage = mode.calculate_voltage(start, end)
    current = mode.calculate_current(elements="cond")
    flux = _flux(mode)

    z_vi = mode.calculate_impedance("VI", start=start, end=end, elements="cond")
    z_pv = mode.calculate_impedance("PV", start=start, end=end)
    z_pi = mode.calculate_impedance("PI", elements="cond")
    np.testing.assert_allclose(z_vi, voltage / current, rtol=1e-12)
    np.testing.assert_allclose(z_pv, abs(voltage) ** 2 / np.conj(flux), rtol=1e-12)
    np.testing.assert_allclose(z_pi, flux / abs(current) ** 2, rtol=1e-12)

    # calculate_power is 0.5 * int (E x H* + E* x H) . z (complex64), i.e. Re(S)
    np.testing.assert_allclose(flux.real, mode.calculate_power().real, rtol=1e-6)

    # Z_VI, Z_PV and Z_PI do not depend on the scaling of the mode
    scale = 2 - 3j
    scaled = dataclasses.replace(mode, E=scale * mode.E, H=scale * mode.H)
    np.testing.assert_allclose(
        scaled.calculate_impedance("VI", start=start, end=end, elements="cond"), z_vi, rtol=1e-10
    )
    np.testing.assert_allclose(
        scaled.calculate_impedance("PV", start=start, end=end), z_pv, rtol=1e-10
    )
    np.testing.assert_allclose(scaled.calculate_impedance("PI", elements="cond"), z_pi, rtol=1e-10)

    _raises(ValueError, mode.calculate_impedance, "ZZ")
    _raises(ValueError, mode.calculate_impedance, "PV")
    _raises(ValueError, mode.calculate_impedance, "VI", start=start, end=end)
    _raises(ValueError, mode.calculate_impedance, "PI")


def _coax_mode(epsilon_r, radius_inner=0.512, radius_outer=2.23039):
    core = shapely.Point(0, 0).buffer(radius_inner)
    isolator = shapely.Point(0, 0).buffer(radius_outer)
    polygons = OrderedDict(
        surface=shapely.LineString(isolator.exterior),
        surface_core=shapely.LinearRing(core.exterior),
        isolator=isolator - core,
    )
    resolutions = dict(isolator={"resolution": 0.2, "distance": 1})
    mesh = from_meshio(
        mesh_from_OrderedDict(
            polygons, resolutions, default_resolution_max=1, mesh_scaling_factor=1e-3
        )
    )
    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex) + epsilon_r
    modes = compute_modes(
        basis0,
        epsilon,
        wavelength=scipy.constants.speed_of_light / 10e9,
        mu_r=1,
        num_modes=1,
        order=2,
        metallic_boundaries=("surface", "surface_core"),
    )
    return modes[0]


def test_coax_tem_consistency():
    radius_inner, radius_outer = 0.512, 2.23039
    start, end = (radius_inner * 1e-3, 0.0), (radius_outer * 1e-3, 0.0)
    for epsilon_r in (1.29, 1.29 * (1 - 0.05j)):
        mode = _coax_mode(epsilon_r, radius_inner, radius_outer)
        analytic = 60 / np.sqrt(epsilon_r) * np.log(radius_outer / radius_inner)

        z_vi = mode.calculate_impedance("VI", start=start, end=end, facets="surface_core")
        z_pv = mode.calculate_impedance("PV", start=start, end=end)
        z_pi = mode.calculate_impedance("PI", facets="surface_core")
        for z in (z_vi, z_pv, z_pi):
            np.testing.assert_allclose(z, analytic, rtol=1e-2)
            # catches a missing conjugate or a real-only power
            atol = 0.05 * abs(analytic.imag) + 1e-3  # analytic.imag is 0 without loss
            np.testing.assert_allclose(z.imag, analytic.imag, rtol=0, atol=atol)
        assert z_vi.real > 0

        # P = V I* / 2
        voltage = mode.calculate_voltage(start, end)
        current = mode.calculate_current(facets="surface_core")
        np.testing.assert_allclose(voltage * np.conj(current) / _flux(mode), 1, rtol=1e-2)

        # the current in the shield is the opposite of the current in the core
        np.testing.assert_allclose(mode.calculate_current(facets="surface"), -current, rtol=1e-2)


def _mesh_microstrip(wsim, hclad, hsi, wcore, hcore):
    # same geometry as docs/electronics/examples/microstrip_waveguide_vary_width.py
    core = box(-wcore / 2, -hcore / 2, wcore / 2, hcore / 2)
    clad = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 + hclad)
    silicon = box(-wsim / 2, -hcore / 2, wsim / 2, -hcore / 2 - hsi)
    combined = shapely.ops.unary_union([clad, silicon])
    polygons = OrderedDict(
        surface=LineString(combined.exterior),
        core_interface=core.buffer(hcore / 10).exterior,
        core=core,
        clad=clad,
        silicon=silicon,
    )
    resolutions = dict(
        core={"resolution": min(wcore / 20, hcore / 4), "distance": 1, "SizeMax": 0.5},
        core_interface={"resolution": min(wcore / 20, hcore / 4), "distance": 1, "SizeMax": 0.5},
    )
    return mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=5)


def test_microstrip_example_regression():
    width, hcore, hsi = 0.4, 0.005, 0.64
    frequency = 10e9
    mesh = from_meshio(_mesh_microstrip(wsim=10, hclad=10, hsi=hsi, wcore=width, hcore=hcore))
    mesh = mesh.scaled((1e-3,) * 2)

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex)
    epsilon[basis0.get_dofs(elements="silicon")] = 9.9 + 0.0005
    epsilon[basis0.get_dofs(elements="clad")] = 1.0
    epsilon[basis0.get_dofs(elements="core")] = 1 - 1j * 1 / (
        18e-6 * 1e-3
    ) / scipy.constants.epsilon_0 / (2 * np.pi * frequency)
    mode = compute_modes(
        basis0,
        epsilon,
        wavelength=scipy.constants.speed_of_light / frequency,
        num_modes=1,
        metallic_boundaries=True,
    )[0]

    # current as computed by the example before it used calculate_impedance
    @Functional(dtype=np.complex64)
    def current_form(w):
        return inner(np.array([w.n[1], -w.n[0]]), w.H)

    (ht, ht_basis), _ = mode.basis.split(mode.H)
    facet_basis = ht_basis.boundary(facets=mesh.boundaries["core_interface"])
    current_example = abs(current_form.assemble(facet_basis, H=facet_basis.interpolate(ht)))

    contour = box(-width / 2, -hcore / 2, width / 2, hcore / 2).buffer(hcore / 10)
    z_pi = mode.calculate_impedance(
        "PI", elements=lambda x: shapely.contains_xy(contour, x[0] * 1e3, x[1] * 1e3)
    )
    np.testing.assert_allclose(z_pi.real, 1 / current_example**2, rtol=1e-5)
    # 61.39 Ohm was measured with the example (mesh dependent)
    np.testing.assert_allclose(z_pi.real, 61.39, rtol=5e-2)
