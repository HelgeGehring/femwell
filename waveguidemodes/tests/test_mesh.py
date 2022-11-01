from typing import Dict, Optional, Tuple, List

import numpy as np
import pygmsh
import gmsh
import shapely
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.ops import split, linemerge

from collections import OrderedDict

from waveguidemodes.mesh import mesh_from_polygons

def geometry(wsim, hclad, hbox, wcore, hcore, offset_core):
    core = Polygon([
            Point(-wcore/2, -hcore/2 + offset_core),
            Point(-wcore/2, hcore/2 + offset_core),
            Point(wcore/2, hcore/2 + offset_core),
            Point(wcore/2, -hcore/2 + offset_core),
        ])
    clad = Polygon([
            Point(-wsim/2, -hcore/2),
            Point(-wsim/2, -hcore/2 + hclad),
            Point(wsim/2, -hcore/2 + hclad),
            Point(wsim/2, -hcore/2),
        ])
    box = Polygon([
            Point(-wsim/2, -hcore/2),
            Point(-wsim/2, -hcore/2 - hbox),
            Point(wsim/2, -hcore/2 - hbox),
            Point(wsim/2, -hcore/2),
        ])

    shapes = OrderedDict()
    shapes["core"] = core 
    shapes["clad"] = clad
    shapes["box"] = box

    return shapes

def test_shared_edge():
    shapes = geometry(wsim = 2, 
                            hclad = 2, 
                            hbox = 2, 
                            wcore = 0.5, 
                            hcore = 0.22, 
                            offset_core = 0
                            )
    mesh = mesh_from_polygons(shapes, resolutions = {})
    assert True

def test_breaking_edge():
    shapes = geometry(wsim = 2, 
                            hclad = 2, 
                            hbox = 2, 
                            wcore = 0.5, 
                            hcore = 0.22, 
                            offset_core = -0.1
                            )
    mesh = mesh_from_polygons(shapes, resolutions = {})
    assert True

def test_inclusion():
    shapes = geometry(wsim = 2, 
                            hclad = 2, 
                            hbox = 2, 
                            wcore = 0.5, 
                            hcore = 0.22, 
                            offset_core = 0.1
                            )
    mesh = mesh_from_polygons(shapes, resolutions = {})
    assert True 

# def test_geometry_with_hole():
#     outer = Polygon([
#             Point(-5, -5),
#             Point(-5, 5),
#             Point(5, 5),
#             Point(5, -5),
#         ], 
#         [Point(-2, -2),
#         Point(-2, 2),
#         Point(2, 2),
#         Point(2, -2),
#         ], 
#         )
#     inner = Polygon([
#             Point(-5, -5),
#             Point(-5, 5),
#             Point(5, 5),
#             Point(5, -5),
#         ])

#     shapes = OrderedDict()
#     shapes["outer"] = outer 
#     shapes["inner"] = inner
#     assert True 