from femwell.polygon_mode_solver import PolygonMaterial, solve_polygon_modes


def test_solve_polygon_modes_returns_sampled_fields():
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
        num_modes=1,
        order=1,
        default_resolution_max=0.5,
    )

    assert len(modes) == 1
    assert modes[0]["n_eff"]["real"] > 0
    assert len(modes[0]["samples"]) == 2
    assert set(modes[0]["samples"][0]["E"]) == {"x", "y", "z"}
    assert set(modes[0]["samples"][0]["H"]["x"]) == {"real", "imag"}


def test_solve_polygon_modes_allows_metallic_boundaries():
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
        sample_points=[(0.0, 0.0)],
        metallic_boundaries=True,
        num_modes=1,
        order=1,
        default_resolution_max=0.5,
    )

    assert modes[0]["n_eff"]["real"] > 0
    assert modes[0]["samples"][0]["point"] == [0.0, 0.0]


def test_solve_polygon_modes_accepts_pec_material_polygon():
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
            ),
            PolygonMaterial(
                name="pec_insert",
                points=[
                    (-0.75, -0.75),
                    (-0.75, -0.55),
                    (-0.55, -0.55),
                    (-0.55, -0.75),
                ],
                epsilon="PEC",
            ),
        ],
        background_epsilon=1.0,
        wavelength=1.55,
        sample_points=[(0.0, 0.0)],
        num_modes=1,
        order=1,
        default_resolution_max=0.5,
    )

    assert modes[0]["n_eff"]["real"] > 0
    assert modes[0]["samples"][0]["point"] == [0.0, 0.0]
