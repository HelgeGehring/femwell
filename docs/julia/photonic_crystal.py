from collections import OrderedDict

import numpy as np
from meshwell.model import Model
from meshwell.prism import Prism
from shapely import Polygon, box
from shapely.geometry import MultiPolygon, Point

if __name__ == "__main__":
    model = Model()

    """
    DEFINE ENTITIES
    """

    a_start = 0.43  # starting periodicity
    a_end = 0.33  # ending periodicity
    s_cav = 0.145  # cavity length
    r = 0.28  # hole radius  (units of a)
    core_thickness = h = 0.22  # waveguide height
    core_width = w = 0.5  # waveguide width

    dair = 1.00  # air padding
    pml_thickness = dpml = 1.00  # PML thickness

    Ndef = 3  # number of defect periods
    a_taper = np.linspace(a_start, a_end, Ndef + 2)
    dgap = a_end - 2 * r * a_end

    Nwvg = 8  # number of waveguide periods
    simulation_length = sx = 2 * (Nwvg * a_start + sum(a_taper)) - dgap + s_cav
    simulation_width = dair + w + dair
    simulation_height = sz = dair + h + dair

    buffer_resolution = 4

    holes = []
    for mm in range(Nwvg):
        holes.append(
            Point((-0.5 * sx + 0.5 * a_start + mm * a_start, 0, 0)).buffer(
                r * a_start, resolution=buffer_resolution
            )
        )
        holes.append(
            Point((+0.5 * sx - 0.5 * a_start - mm * a_start, 0, 0)).buffer(
                r * a_start, resolution=buffer_resolution
            )
        )

    for mm in range(Ndef + 2):
        holes.append(
            Point(
                (
                    -0.5 * sx
                    + Nwvg * a_start
                    + (sum(a_taper[:mm]) if mm > 0 else 0)
                    + 0.5 * a_taper[mm],
                    0,
                    0,
                )
            ).buffer(r * a_taper[mm], resolution=buffer_resolution)
        )
        holes.append(
            Point(
                (
                    +0.5 * sx
                    - Nwvg * a_start
                    - (sum(a_taper[:mm]) if mm > 0 else 0)
                    - 0.5 * a_taper[mm],
                    0,
                    0,
                )
            ).buffer(r * a_taper[mm], resolution=buffer_resolution)
        )

    # Create a multipolygon from the holes
    holes_multipolygon = MultiPolygon(holes)

    # Core
    core_polygon = Polygon(
        shell=(
            (-simulation_length / 2 - dpml, -core_width / 2),
            (simulation_length / 2 + dpml, -core_width / 2),
            (simulation_length / 2 + dpml, core_width / 2),
            (-simulation_length / 2 - dpml, core_width / 2),
        ),
    )
    core_buffers = {
        0: 0.0,
        core_thickness: 0.0,
    }
    core = Prism(
        polygons=core_polygon - holes_multipolygon,
        buffers=core_buffers,
        model=model,
    )

    cladding_polygon = box(
        xmin=-simulation_length / 2,
        ymin=-simulation_width / 2,
        xmax=simulation_length / 2,
        ymax=simulation_width / 2,
    )

    pml_polygon = box(
        xmin=-simulation_length / 2 - pml_thickness,
        ymin=-simulation_width / 2 - pml_thickness,
        xmax=simulation_length / 2 + pml_thickness,
        ymax=simulation_width / 2 + pml_thickness,
    )

    PML = Prism(
        polygons=pml_polygon,
        buffers={
            -simulation_height - pml_thickness: 0.0,
            simulation_height + pml_thickness: 0.0,
        },
        model=model,
    )

    """
    ASSEMBLE AND NAME ENTITIES
    """
    entities = OrderedDict()
    entities["core"] = core
    entities["PML"] = PML

    mesh = model.mesh(
        entities_dict=entities, verbosity=0, filename="mesh.msh", default_characteristic_length=1
    )
