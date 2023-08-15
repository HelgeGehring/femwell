from collections import OrderedDict

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
    device_polygon = Polygon(shell=((0.25, 0.25), (0.25, 0.75), (0.75, 0.75), (0.25, 0.25)))
    device_buffers = {
        0: 0.0,
        device_thickness: 0.0,
    }

    device = Prism(
        polygons=device_polygon,
        buffers={
            0: 0.0,
            device_thickness: 0.0,
        },
        model=model,
    )

    background_polygon = box(xmin=0, ymin=0, xmax=simulation_length, ymax=simulation_width)

    substrate = Prism(
        polygons=background_polygon,
        buffers={
            -substrate_thickness: 0.0,
            0: 0.0,
        },
        model=model,
    )

    cladding = Prism(
        polygons=background_polygon,
        buffers={
            0: 0.0,
            cladding_thickness: 0.0,
        },
        model=model,
    )

    PML_cladding = Prism(
        polygons=background_polygon,
        buffers={
            cladding_thickness: 0.0,
            cladding_thickness + pml_thickness: 0.0,
        },
        model=model,
    )

    PML_substrate = Prism(
        polygons=background_polygon,
        buffers={
            -substrate_thickness: 0.0,
            -substrate_thickness - pml_thickness: 0.0,
        },
        model=model,
    )

    """
    BOUNDARIES
    """

    right_polygon = box(
        xmin=simulation_length,
        ymin=0,
        xmax=simulation_length + 1,
        ymax=simulation_width,
    )
    right = Prism(
        polygons=right_polygon,
        buffers={
            -substrate_thickness - pml_thickness: 0.0,
            cladding_thickness + pml_thickness: 0.0,
        },
        model=model,
    )

    left_polygon = box(xmin=-1, ymin=0, xmax=0, ymax=simulation_width)
    left = Prism(
        polygons=left_polygon,
        buffers={
            -substrate_thickness - pml_thickness: 0.0,
            cladding_thickness + pml_thickness: 0.0,
        },
        model=model,
    )

    up_polygon = box(
        xmin=0, ymin=simulation_width, xmax=simulation_length, ymax=simulation_width + 1
    )
    up = Prism(
        polygons=up_polygon,
        buffers={
            -substrate_thickness - pml_thickness: 0.0,
            cladding_thickness + pml_thickness: 0.0,
        },
        model=model,
    )

    down_polygon = box(xmin=0, ymin=-1, xmax=simulation_length, ymax=0)
    down = Prism(
        polygons=down_polygon,
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
    entities["substrate"] = substrate
    # entities["device"] = device
    entities["cladding"] = cladding
    entities["PML_cladding"] = PML_cladding
    entities["PML_substrate"] = PML_substrate

    boundary_entities = {}
    boundary_entities["right"] = right
    boundary_entities["left"] = left
    boundary_entities["up"] = up
    boundary_entities["down"] = down

    mesh = model.mesh(
        entities_dict=entities,
        boundaries_dict=boundary_entities,
        verbosity=0,
        filename="mesh.msh",
        periodic_entities=[
            (x + "___" + s1, x + "___" + s2)
            for x in ("cladding", "substrate", "PML_cladding", "PML_substrate")
            for (s1, s2) in (("left", "right"), ("up", "down"))
        ],
        default_characteristic_length=0.2,
        resolutions={"substrate": {"resolution": 0.07}, "cladding": {"resolution": 0.07}},
    )
