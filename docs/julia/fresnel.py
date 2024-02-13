from shapely import Polygon, box
from meshwell.prism import Prism
from meshwell.model import Model

if __name__ == "__main__":
    model = Model()

    """
    DEFINE ENTITIES
    """

    simvolx = 60
    simvoly = 60
    simvolz = 30
    pml_z = 1

    vol1_polygon = box(xmin = -simvolx/2, ymin = -simvoly/2, xmax = simvolx/2, ymax = simvoly/2)

    vol1 = Prism(
        polygons=vol1_polygon,
        buffers={
            -simvolz/2: 0.0,
            simvolz/2: 0.0,
        },
        model=model,
        physical_name="volume1",
        mesh_order=3,
        resolution={"resolution": 0.5, "SizeMax": 1.0, "DistMax": 1.0},
    )

    PML_bot = Prism(
        polygons=vol1_polygon,
        buffers={
            -simvolz/2 - pml_z: 0.0,
            -simvolz/2: 0.0,
        },
        model=model,
        physical_name="PML_bottom",
        mesh_order=2,
    )

    PML_top = Prism(
        polygons=vol1_polygon,
        buffers={
            simvolz/2: 0.0,
            simvolz/2 + pml_z: 0.0,
        },
        model=model,
        physical_name="PML_top",
        mesh_order=1,
    )

    """
    BOUNDARIES
    """

    right_polygon = box(xmin = simvolx/2, ymin = -simvoly/2, xmax = simvolx/2 + 1, ymax = simvoly/2)
    left_polygon = box(xmin = -simvolx/2 - 1, ymin = -simvoly/2, xmax = -simvolx/2, ymax = simvoly/2)
    front_polygon = box(xmin = -simvolx/2, ymin = simvoly/2, xmax = simvolx/2, ymax = simvoly/2 + 1)
    back_polygon = box(xmin = -simvolx/2, ymin = -simvoly/2 - 1, xmax = simvolx/2, ymax = -simvoly/2)

    boundary_buffers = {-simvolz/2 - pml_z: 0.0, simvolz/2 + pml_z: 0.0}

    right = Prism(
        polygons=right_polygon,
        buffers=boundary_buffers,
        model=model,
        physical_name="right",
        mesh_bool=False,
        mesh_order=0
    )

    left = Prism(
        polygons=left_polygon,
        buffers=boundary_buffers,
        model=model,
        physical_name="left",
        mesh_bool=False,
        mesh_order=0
    )

    up = Prism(
        polygons=front_polygon,
        buffers=boundary_buffers,
        model=model,
        physical_name="up",
        mesh_bool=False,
        mesh_order=0
    )

    down = Prism(
        polygons=back_polygon,
        buffers=boundary_buffers,
        model=model,
        physical_name="down",
        mesh_bool=False,
        mesh_order=0
    )

    """
    ASSEMBLE AND NAME ENTITIES
    """
    entities = [
        vol1,
        PML_bot,
        PML_top,
        up,
        down,
        left,
        right
    ]

    mesh = model.mesh(
        entities_list=entities,
        verbosity=0,
        global_scaling = 1e-6,
        filename="mesh.msh",
        periodic_entities=[
            (x + "___" + s1, x + "___" + s2)
            for x in ("volume1", "PML_bot", "PML_top")
            for (s1, s2) in (("left", "right"), ("up", "down"))
        ]
    )
