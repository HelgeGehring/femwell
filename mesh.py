import pygmsh
from pygmsh.common.polygon import Polygon

um = 1

with pygmsh.geo.Geometry() as geom:
    poly = Polygon(geom,
                   [
                       [-0.25*um, -0.125*um],
                       [0.25*um, -0.125*um],
                       [0.25*um, 0.125*um],
                       [-0.25*um, 0.125*um],
                   ],
                   mesh_size=0.02*um
                   )
    geom.add_physical([poly], 'Core')
    # geom.add_physical(poly.curves, 'Interface')

    poly2 = Polygon(geom,
                    [
                        [-1.2*um, -.6*um],
                        [1.2*um, -.6*um],
                        [1.2*um, .6*um],
                        [-1.2*um, .6*um],
                    ],
                    mesh_size=0.05*um,
                    holes=[poly]
                    )
    geom.add_physical([poly2], 'Cladding')

    # field0 = geom.add_boundary_layer(
    #    edges_list=[poly.curves[0]],
    #   lcmin=0.02*um,
    #    lcmax=0.2*um,
    #    distmin=0.02*um,
    #    distmax=0.2*um,
    # )
    # field1 = geom.add_boundary_layer(
    #     nodes_list=[poly.points[2]],
    #     lcmin=0.05,
    #     lcmax=0.2,
    #     distmin=0.1,
    #     distmax=0.4,
    # )
    # geom.set_background_mesh([field0, ], operator="Min")
    mesh = geom.generate_mesh(dim=2)
    print(mesh.cells)
    print(mesh.point_data)

    import gmsh

    gmsh.write('mesh.msh')
