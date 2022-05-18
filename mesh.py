import pygmsh
from pygmsh.common.polygon import Polygon

with pygmsh.geo.Geometry() as geom:
    poly = Polygon(geom,
                   [
                       [-0.25e-6, -0.11e-6],
                       [0.25e-6, -0.11e-6],
                       [0.25e-6, 0.11e-6],
                       [-0.25e-6, 0.11e-6],
                   ],
                   mesh_size=0.04e-6
                   )
    geom.add_physical([poly], 'Core')
    # geom.add_physical(poly.curves, 'Interface')

    poly2 = Polygon(geom,
                    [
                        [-0.6e-6, -0.3e-6],
                        [0.6e-6, -0.3e-6],
                        [0.6e-6, .3e-6],
                        [-0.6e-6, .3e-6],
                    ],
                    mesh_size=0.04e-6,
                    holes=[poly]
                    )
    geom.add_physical([poly2], 'Cladding')

    # field0 = geom.add_boundary_layer(
    #    edges_list=[poly.curves[0]],
    #   lcmin=0.02e-6,
    #    lcmax=0.2e-6,
    #    distmin=0.02e-6,
    #    distmax=0.2e-6,
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
