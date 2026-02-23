from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import gmsh
import numpy as np
import pygmsh
import shapely
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.ops import linemerge, split

from femwell.mesh import mesh_from_OrderedDict


def geometry(wsim, hclad, hbox, wcore, hcore, offset_core):
    core = Polygon(
        [
            (-wcore / 2, -hcore / 2 + offset_core),
            (-wcore / 2, hcore / 2 + offset_core),
            (wcore / 2, hcore / 2 + offset_core),
            (wcore / 2, -hcore / 2 + offset_core),
        ]
    )
    clad = Polygon(
        [
            (-wsim / 2, -hcore / 2),
            (-wsim / 2, -hcore / 2 + hclad),
            (wsim / 2, -hcore / 2 + hclad),
            (wsim / 2, -hcore / 2),
        ]
    )
    box = Polygon(
        [
            (-wsim / 2, -hcore / 2),
            (-wsim / 2, -hcore / 2 - hbox),
            (wsim / 2, -hcore / 2 - hbox),
            (wsim / 2, -hcore / 2),
        ]
    )

    shapes = OrderedDict()
    shapes["core"] = core
    shapes["clad"] = clad
    shapes["box"] = box

    return shapes


def test_shared_edge():
    shapes = geometry(wsim=2, hclad=2, hbox=2, wcore=0.5, hcore=0.22, offset_core=0)
    mesh = mesh_from_OrderedDict(shapes, resolutions={})


def test_breaking_edge():
    shapes = geometry(wsim=2, hclad=2, hbox=2, wcore=0.5, hcore=0.22, offset_core=-0.1)
    mesh = mesh_from_OrderedDict(shapes, resolutions={})


def test_inclusion():
    shapes = geometry(wsim=2, hclad=2, hbox=2, wcore=0.5, hcore=0.22, offset_core=0.1)
    mesh = mesh_from_OrderedDict(shapes, resolutions={})


def test_lines():
    wmode = 1
    wsim = 2
    hclad = 2
    hbox = 2
    wcore = 0.5
    hcore = 0.22
    offset_core = -0.1
    offset_core2 = 1

    # Lines can be added, which is useful to define boundary conditions at various simulation edges
    left_edge = LineString([(-wsim / 2, -hcore / 2 - hbox), (-wsim / 2, -hcore / 2 + hclad)])
    right_edge = LineString([(wsim / 2, -hcore / 2 - hbox), (wsim / 2, -hcore / 2 + hclad)])
    top_edge = LineString([(-wsim / 2, -hcore / 2 + hclad), (wsim / 2, -hcore / 2 + hclad)])
    bottom_edge = LineString([(-wsim / 2, -hcore / 2 - hbox), (wsim / 2, -hcore / 2 - hbox)])

    # The order in which objects are inserted into the OrderedDict determines overrrides
    shapes = OrderedDict()
    shapes["left_edge"] = left_edge
    shapes["right_edge"] = right_edge
    shapes["top_edge"] = top_edge
    shapes["bottom_edge"] = bottom_edge

    shapes.update(
        geometry(
            wsim=wsim,
            hclad=hclad,
            hbox=hbox,
            wcore=wcore,
            hcore=hcore,
            offset_core=offset_core,
        )
    )
    mesh = mesh_from_OrderedDict(shapes, resolutions={})
