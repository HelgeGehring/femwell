from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

import gmsh
import numpy as np
import pygmsh
import shapely
from shapely.geometry import (
    LinearRing,
    LineString,
    MultiLineString,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import linemerge, split


class MeshTracker:
    def __init__(self, model, atol=1e-3):
        """
        Map between shapely and gmsh
        Shapely is useful for built-in geometry equivalencies and extracting orientation, instead of doing it manually
        We can also track information about the entities using labels (useful for selective mesh refinement later)
        """
        self.shapely_points = []
        self.gmsh_points = []
        self.points_labels = []
        self.shapely_xy_segments = []
        self.gmsh_xy_segments = []
        self.xy_segments_main_labels = []
        self.xy_segments_secondary_labels = []
        self.shapely_xy_surfaces = []
        self.gmsh_xy_surfaces = []
        self.xy_surfaces_labels = []
        self.model = model
        self.atol = atol

    """
    Retrieve existing geometry
    """

    def get_point_index(self, xy_point):
        for index, shapely_point in enumerate(self.shapely_points):
            if xy_point.equals_exact(shapely_point, self.atol):
                return index
        return None

    def get_xy_segment_index_and_orientation(self, xy_point1, xy_point2):
        xy_line = shapely.geometry.LineString([xy_point1, xy_point2])
        for index, shapely_line in enumerate(self.shapely_xy_segments):
            if xy_line.equals(shapely_line):
                first_xy_line, last_xy_line = xy_line.boundary.geoms
                first_xy, last_xy = shapely_line.boundary.geoms
                return (index, first_xy_line.equals(first_xy))
        return None, 1

    def get_gmsh_points_from_label(self, label):
        indices = [idx for idx, value in enumerate(self.points_labels) if value == label]
        return [self.gmsh_points[index]._id for index in indices]

    def get_gmsh_xy_lines_from_label(self, label):
        indices = [idx for idx, value in enumerate(self.xy_segments_main_labels) if value == label]
        return [self.gmsh_xy_segments[index]._id for index in indices]

    def get_gmsh_xy_surfaces_from_label(self, label):
        indices = [idx for idx, value in enumerate(self.xy_surfaces_labels) if value == label]
        return [self.gmsh_xy_surfaces[index]._id for index in indices]

    """
    Channel loop utilities (no need to track)
    """

    def xy_channel_loop_from_vertices(self, vertices, label):
        edges = []
        for vertex1, vertex2 in [(vertices[i], vertices[i + 1]) for i in range(len(vertices) - 1)]:
            gmsh_line, orientation = self.add_get_xy_segment(vertex1, vertex2, label)
            if orientation:
                edges.append(gmsh_line)
            else:
                edges.append(-gmsh_line)
        return self.model.add_curve_loop(edges)

    """
    Adding geometry
    """

    def add_get_point(self, shapely_xy_point, label=None):
        """
        Add a shapely point to the gmsh model, or retrieve the existing gmsh model points with equivalent coordinates (within tol.)

        Args:
            shapely_xy_point (shapely.geometry.Point): x, y coordinates
            resolution (float): gmsh resolution at that point
        """
        index = self.get_point_index(shapely_xy_point)
        if index is not None:
            gmsh_point = self.gmsh_points[index]
        else:
            gmsh_point = self.model.add_point([shapely_xy_point.x, shapely_xy_point.y])
            self.shapely_points.append(shapely_xy_point)
            self.gmsh_points.append(gmsh_point)
            self.points_labels.append(label)
        return gmsh_point

    def add_get_xy_segment(self, shapely_xy_point1, shapely_xy_point2, label):
        """
        Add a shapely segment (2-point line) to the gmsh model in the xy plane, or retrieve the existing gmsh segment with equivalent coordinates (within tol.)

        Args:
            shapely_xy_point1 (shapely.geometry.Point): first x, y coordinates
            shapely_xy_point2 (shapely.geometry.Point): second x, y coordinates
        """
        index, orientation = self.get_xy_segment_index_and_orientation(
            shapely_xy_point1, shapely_xy_point2
        )
        if index is not None:
            gmsh_segment = self.gmsh_xy_segments[index]
            self.xy_segments_secondary_labels[index] = label
        else:
            gmsh_segment = self.model.add_line(
                self.add_get_point(shapely_xy_point1),
                self.add_get_point(shapely_xy_point2),
            )
            self.shapely_xy_segments.append(
                shapely.geometry.LineString([shapely_xy_point1, shapely_xy_point2])
            )
            self.gmsh_xy_segments.append(gmsh_segment)
            self.xy_segments_main_labels.append(label)
            self.xy_segments_secondary_labels.append(None)
        return gmsh_segment, orientation

    def add_get_xy_line(self, shapely_xy_curve, label):
        """
        Add a shapely line (multi-point line) to the gmsh model in the xy plane, or retrieve the existing gmsh segment with equivalent coordinates (within tol.)

        Args:
            shapely_xy_curve (shapely.geometry.LineString): curve
        """
        segments = []
        for shapely_xy_point1, shapely_xy_point2 in zip(
            shapely_xy_curve.coords[:-1], shapely_xy_curve.coords[1:]
        ):
            gmsh_segment, orientation = self.add_get_xy_segment(
                Point(shapely_xy_point1), Point(shapely_xy_point2), label
            )
            if orientation:
                segments.append(gmsh_segment)
            else:
                segments.append(-gmsh_segment)
        self.model.add_physical(segments, f"{label}")

    def add_xy_surface(self, shapely_xy_polygons, label=None, physical=True):
        """
        Add a xy surface corresponding to shapely_xy_polygon, or retrieve the existing gmsh model surface with equivalent coordinates (within tol.)

        Args:
            shapely_xy_polygons (shapely.geometry.(Multi)Polygon):
        """
        # Create surface
        surfaces_to_label = []

        for shapely_xy_polygon in (
            shapely_xy_polygons.geoms
            if hasattr(shapely_xy_polygons, "geoms")
            else [shapely_xy_polygons]
        ):
            hole_loops = []
            exterior_vertices = []
            # Parse holes
            for polygon_hole in list(shapely_xy_polygon.interiors):
                hole_vertices = []
                for vertex in shapely.geometry.MultiPoint(polygon_hole.coords).geoms:
                    # gmsh_point = self.add_get_point(vertex, label)
                    hole_vertices.append(vertex)
                hole_loops.append(self.xy_channel_loop_from_vertices(hole_vertices, label))
            # Parse boundary
            for vertex in shapely.geometry.MultiPoint(shapely_xy_polygon.exterior.coords).geoms:
                # gmsh_point = self.add_get_point(vertex, label)
                exterior_vertices.append(vertex)
            channel_loop = self.xy_channel_loop_from_vertices(exterior_vertices, label)

            # Create and log surface
            gmsh_surface = self.model.add_plane_surface(channel_loop, holes=hole_loops)
            self.shapely_xy_surfaces.append(shapely_xy_polygon)
            self.gmsh_xy_surfaces.append(gmsh_surface)
            self.xy_surfaces_labels.append(label)
            surfaces_to_label.append(gmsh_surface)
        if physical:
            self.model.add_physical(surfaces_to_label, f"{label}")
