from collections import OrderedDict
from typing import Dict
from webbrowser import BackgroundBrowser

from meshwell.gmsh_entity import GMSH_entity
from meshwell.model import Model
from meshwell.prism import Prism
from shapely import Polygon, box

if __name__ == "__main__":
    model = Model()

    """
    DEFINE ENTITIES
    """

    simulation_width = 1
    simulation_length = 1

    substrate_thickness = 1
    cladding_thickness = 1
    pml_thickness = 1

    device_thickness = 0.1
    device_polygon = Polygon(shell=((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.0, 0.0)))
    device = Prism(
        polygons=device_polygon,
        buffers={0: 0.0, device_thickness: 0.0},
        model=model,
    )

    background_polygon = box(xmin=0, ymin=0, xmax=simulation_length, ymax=simulation_width)

    substrate = Prism(
        polygons=background_polygon,
        buffers={-substrate_thickness: 0.0, 0: 0.0},
        model=model,
    )

    cladding = Prism(
        polygons=background_polygon,
        buffers={0: 0.0, cladding_thickness: 0.0},
        model=model,
    )

    PML_cladding = Prism(
        polygons=background_polygon,
        buffers={cladding_thickness: 0.0, cladding_thickness + pml_thickness: 0.0},
        model=model,
    )

    PML_substrate = Prism(
        polygons=background_polygon,
        buffers={-substrate_thickness: 0.0, -substrate_thickness - pml_thickness: 0.0},
        model=model,
    )
    import shapely.ops

    boundary = Prism(
        polygons=shapely.ops.unary_union([box(1, 0, 2, 1), box(-1, 0, 0, 1)]),
        buffers={
            -substrate_thickness - pml_thickness: 0.0,
            cladding_thickness + pml_thickness: 0.0,
        },
        model=model,
    )

    """
    ASSEMBLE AND NAME ENTITIES
    """
    entities = OrderedDict()
    # entities["device"] = device
    entities["substrate"] = substrate
    entities["cladding"] = cladding
    entities["PML_cladding"] = PML_cladding
    entities["PML_substrate"] = PML_substrate

    boundaries = dict(boundary=boundary)

    mesh = model.mesh(
        entities,
        boundaries_dict=boundaries,
        verbosity=5,
        filename="mesh.msh",
        default_characteristic_length=0.2,
        resolutions={"substrate": {"resolution": 0.07}, "cladding": {"resolution": 0.07}},
    )
