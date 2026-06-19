"""Solve modes from simple polygon input and print sampled fields as JSON."""

from __future__ import annotations

import json

from femwell.polygon_mode_solver import PolygonMaterial, solve_polygon_modes


def main() -> None:
    modes = solve_polygon_modes(
        bbox=(-1.0, -1.0, 1.0, 1.0),
        materials=[
            PolygonMaterial(
                name="core",
                points=[
                    (-0.2, -0.1),
                    (-0.2, 0.1),
                    (0.2, 0.1),
                    (0.2, -0.1),
                ],
                epsilon=2.25,
                resolution=0.2,
            )
        ],
        background_epsilon=1.0,
        wavelength=1.55,
        sample_points=[(0.0, 0.0), (0.5, 0.5)],
        metallic_boundaries=True,
        num_modes=1,
        order=1,
        default_resolution_max=0.5,
    )
    print(json.dumps(modes, indent=2))


if __name__ == "__main__":
    main()
