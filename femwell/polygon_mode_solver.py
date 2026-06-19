"""Small adapter for polygon-described waveguide mode solves."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np
import shapely
from shapely.geometry import LineString, Polygon
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio

from femwell.maxwell.waveguide import compute_modes
from femwell.mesh import mesh_from_OrderedDict

Point = tuple[float, float]
BBox = tuple[float, float, float, float]


@dataclass(frozen=True)
class PolygonMaterial:
    """A named polygon and relative permittivity for a mode-solver input."""

    name: str
    points: Sequence[Point]
    epsilon: complex
    resolution: float | None = None
    refinement_distance: float = 0.5


@dataclass(frozen=True)
class PECLine:
    """A named perfect-electric-conductor line constraint."""

    name: str
    points: Sequence[Point]


def solve_polygon_modes(
    *,
    bbox: BBox,
    materials: Sequence[PolygonMaterial],
    wavelength: float,
    sample_points: Sequence[Point],
    background_epsilon: complex = 1.0,
    background_name: str = "background",
    pec_lines: Sequence[PECLine] = (),
    metallic_boundaries: bool | Sequence[str] = False,
    num_modes: int = 1,
    order: int = 1,
    default_resolution_max: float = 0.5,
    solver: str = "scipy",
) -> list[dict]:
    """Solve modes from bbox/polygon inputs and return sampled fields.

    The input shape mirrors simple integration payloads: a bounding box, a
    list of named polygons with relative permittivities, optional PEC line
    labels, and sample coordinates where the solved electric and magnetic
    fields are needed.

    Args:
        bbox: ``(xmin, ymin, xmax, ymax)`` simulation bounds.
        materials: Ordered material polygons. Earlier polygons override later
            background material when the mesh is tiled.
        wavelength: Vacuum wavelength in the same length unit as the geometry.
        sample_points: Coordinates where E/H fields are returned.
        background_epsilon: Relative permittivity for the bbox remainder.
        background_name: Mesh domain name for the bbox remainder.
        pec_lines: Optional line labels to use as metallic boundaries.
        metallic_boundaries: ``True`` for all exterior boundaries, ``False``
            for none, or a list of existing line labels.
        num_modes: Number of eigenmodes to compute.
        order: Finite-element order passed to :func:`compute_modes`.
        default_resolution_max: Maximum gmsh element size.
        solver: Eigenvalue solver passed to :func:`compute_modes`.

    Returns:
        A JSON-serializable list with mode metadata and sampled E/H components.
    """

    shapes, resolutions = _build_shapes(
        bbox=bbox,
        materials=materials,
        background_name=background_name,
        pec_lines=pec_lines,
    )

    mesh = from_meshio(
        mesh_from_OrderedDict(
            shapes,
            resolutions,
            default_resolution_max=default_resolution_max,
        )
    )
    basis_epsilon = Basis(mesh, ElementTriP0())
    epsilon = _build_epsilon_field(
        basis_epsilon=basis_epsilon,
        materials=materials,
        background_name=background_name,
        background_epsilon=background_epsilon,
    )

    metallic = _metallic_boundary_argument(metallic_boundaries, pec_lines)
    modes = compute_modes(
        basis_epsilon,
        epsilon,
        wavelength=wavelength,
        metallic_boundaries=metallic,
        num_modes=num_modes,
        order=order,
        solver=solver,
    )

    coordinates = _coordinates_array(sample_points)
    return [_serialize_mode(mode, sample_points, coordinates) for mode in modes]


def _build_shapes(
    *,
    bbox: BBox,
    materials: Sequence[PolygonMaterial],
    background_name: str,
    pec_lines: Sequence[PECLine],
) -> tuple[OrderedDict[str, LineString | Polygon], dict[str, dict[str, float]]]:
    shapes: OrderedDict[str, LineString | Polygon] = OrderedDict()
    resolutions: dict[str, dict[str, float]] = {}

    for pec_line in pec_lines:
        shapes[pec_line.name] = LineString(pec_line.points)

    for material in materials:
        shapes[material.name] = Polygon(material.points)
        if material.resolution is not None:
            resolutions[material.name] = {
                "resolution": material.resolution,
                "distance": material.refinement_distance,
            }

    shapes[background_name] = shapely.box(*bbox)
    return shapes, resolutions


def _build_epsilon_field(
    *,
    basis_epsilon: Basis,
    materials: Sequence[PolygonMaterial],
    background_name: str,
    background_epsilon: complex,
) -> np.ndarray:
    epsilon_values = [background_epsilon, *(m.epsilon for m in materials)]
    dtype = complex if _has_complex_values(epsilon_values) else float
    epsilon = basis_epsilon.zeros(dtype=dtype)
    epsilon[basis_epsilon.get_dofs(elements=background_name)] = background_epsilon

    for material in materials:
        epsilon[basis_epsilon.get_dofs(elements=material.name)] = material.epsilon

    return epsilon


def _has_complex_values(values: Iterable[complex]) -> bool:
    return any(np.iscomplexobj(value) and np.imag(value) != 0 for value in values)


def _metallic_boundary_argument(
    metallic_boundaries: bool | Sequence[str],
    pec_lines: Sequence[PECLine],
) -> bool | list[str]:
    if metallic_boundaries is True:
        return True

    boundary_names = list(metallic_boundaries or [])
    boundary_names.extend(pec_line.name for pec_line in pec_lines)
    return boundary_names if boundary_names else False


def _coordinates_array(points: Sequence[Point]) -> np.ndarray:
    return np.asarray(points, dtype=float).T


def _serialize_mode(mode, points: Sequence[Point], coordinates: np.ndarray) -> dict:
    e_values = _sample_field(mode.basis, mode.E, coordinates)
    h_values = _sample_field(mode.basis, mode.H, coordinates)

    return {
        "n_eff": _serialize_complex(mode.n_eff),
        "k": _serialize_complex(mode.k),
        "frequency": float(mode.frequency),
        "samples": [
            {
                "point": [float(point[0]), float(point[1])],
                "E": _serialize_vector(e_values, index),
                "H": _serialize_vector(h_values, index),
            }
            for index, point in enumerate(points)
        ],
    }


def _sample_field(
    basis: Basis,
    field: np.ndarray,
    coordinates: np.ndarray,
) -> Mapping[str, np.ndarray]:
    (field_t, field_t_basis), (field_z, field_z_basis) = basis.split(field)
    values_t = field_t_basis.interpolator(field_t)(coordinates)
    values_z = field_z_basis.interpolator(field_z)(coordinates)

    return {
        "x": values_t[0],
        "y": values_t[1],
        "z": values_z,
    }


def _serialize_vector(values: Mapping[str, np.ndarray], index: int) -> dict:
    return {
        "x": _serialize_complex(values["x"][index]),
        "y": _serialize_complex(values["y"][index]),
        "z": _serialize_complex(values["z"][index]),
    }


def _serialize_complex(value: complex) -> dict[str, float]:
    return {"real": float(np.real(value)), "imag": float(np.imag(value))}
