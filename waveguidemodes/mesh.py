from typing import Dict, Optional, Tuple, List

import numpy as np
import pygmsh
import gmsh
import shapely
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.ops import split, linemerge

from collections import OrderedDict

class MeshTracker():
    def __init__(self, model, atol=1E-3):
        '''
        Map between shapely and gmsh 
        Shapely is useful for built-in geometry equivalencies and extracting orientation, instead of doing it manually
        We can also track information about the entities using labels (useful for selective mesh refinement later)
        '''
        self.shapely_points = []
        self.gmsh_points = []
        self.points_labels = []
        self.shapely_xy_lines = []
        self.gmsh_xy_lines = []
        self.xy_lines_labels = []
        self.gmsh_xy_surfaces = []
        self.xy_surfaces_labels = []
        self.model = model
        self.atol = atol

    """
    Retrieve existing geometry
    """
    def get_point_index(self, xy_point):
        for index, shapely_point in enumerate(self.shapely_points):
            if xy_point.equals_exact(shapely_point, self.atol) :
                return index
        return None

    def get_xy_line_index_and_orientation(self, xy_point1, xy_point2):
        xy_line = shapely.geometry.LineString([xy_point1, xy_point2])
        for index, shapely_line in enumerate(self.shapely_xy_lines):
            if xy_line.equals(shapely_line):
                first_xy_line, last_xy_line = xy_line.boundary.geoms
                first_xy, last_xy = shapely_line.boundary.geoms
                if first_xy_line.equals(first_xy):
                    return index, True
                else:
                    return index, False
        return None, 1

    def get_gmsh_points_from_label(self, label):
        indices = [idx for idx, value in enumerate(self.points_labels) if value == label]
        entities = []
        for index in indices:
            entities.append(self.gmsh_points[index]._id)
        return entities

    def get_gmsh_xy_lines_from_label(self, label):
        indices = [idx for idx, value in enumerate(self.xy_lines_labels) if value == label]
        entities = []
        for index in indices:
            entities.append(self.gmsh_xy_lines[index]._id)
        return entities

    def get_gmsh_xy_surfaces_from_label(self, label):
        indices = [idx for idx, value in enumerate(self.xy_surfaces_labels) if value == label]
        entities = []
        for index in indices:
            entities.append(self.gmsh_xy_surfaces[index]._id)
        return entities

    """
    Channel loop utilities (no need to track)
    """
    def xy_channel_loop_from_vertices(self, vertices, label):
        edges = []
        for vertex1, vertex2 in [(vertices[i], vertices[i + 1]) for i in range(0, len(vertices)-1)]:
            gmsh_line, orientation = self.add_get_xy_line(vertex1, vertex2, label)
            if orientation:
                edges.append(gmsh_line)
            else:
                edges.append(-gmsh_line)
        channel_loop = self.model.add_curve_loop(edges)
        return channel_loop

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


    def add_get_xy_line(self, shapely_xy_point1, shapely_xy_point2, label):
        """
        Add a shapely line to the gmsh model in the xy plane, or retrieve the existing gmsh model line with equivalent coordinates (within tol.)

        Args:
            shapely_xy_point1 (shapely.geometry.Point): first x, y coordinates
            shapely_xy_point2 (shapely.geometry.Point): second x, y coordinates
        """
        index, orientation = self.get_xy_line_index_and_orientation(shapely_xy_point1, shapely_xy_point2)
        if index is not None:
            gmsh_line = self.gmsh_xy_lines[index]
        else:
            gmsh_line = self.model.add_line(self.add_get_point(shapely_xy_point1), self.add_get_point(shapely_xy_point2))
            self.shapely_xy_lines.append(shapely.geometry.LineString([shapely_xy_point1, shapely_xy_point2]))
            self.gmsh_xy_lines.append(gmsh_line)
            self.xy_lines_labels.append(label)
        return gmsh_line, orientation

    def add_xy_surface(self, shapely_xy_polygon, label=None):
        """
        Add a xy surface corresponding to shapely_xy_polygon, or retrieve the existing gmsh model surface with equivalent coordinates (within tol.)

        Args:
            shapely_xy_polygon (shapely.geometry.Polygon):
        """
        # Create surface
        exterior_vertices = []
        hole_loops = []

        # Parse holes
        for polygon_hole in list(shapely_xy_polygon.interiors):
            hole_vertices = []
            for vertex in shapely.geometry.MultiPoint(polygon_hole.coords).geoms:
                gmsh_point = self.add_get_point(vertex, label)
                hole_vertices.append(vertex)
            hole_loops.append(self.xy_channel_loop_from_vertices(hole_vertices, label))
        # Parse boundary
        for vertex in shapely.geometry.MultiPoint(shapely_xy_polygon.exterior.coords).geoms:
            gmsh_point = self.add_get_point(vertex, label)
            exterior_vertices.append(vertex)
        channel_loop = self.xy_channel_loop_from_vertices(exterior_vertices, label)

        # Create and log surface
        gmsh_surface = self.model.add_plane_surface(channel_loop, holes=hole_loops)
        self.gmsh_xy_surfaces.append(gmsh_surface)
        self.xy_surfaces_labels.append(label)

        return gmsh_surface


def mesh_from_polygons(
    polygon_dict: OrderedDict,
    resolutions: Optional[Dict[str, float]] = None,
    default_resolution_min: float = 0.01,
    default_resolution_max: float = 0.1,
    filename: Optional[str] = None,
):

    import gmsh


    with pygmsh.occ.geometry.Geometry() as geometry:

        gmsh.initialize()

        # geometry = pygmsh.occ.geometry.Geometry()
        geometry.characteristic_length_min = default_resolution_min
        geometry.characteristic_length_max = default_resolution_max

        model = geometry.__enter__()

        # Break up surfaces in order so that plane is tiled with non-overlapping layers
        polygons_tiled_dict = OrderedDict()
        for lower_index, (lower_name, lower_polygon) in reversed(list(enumerate(polygon_dict.items()))):
            diff_polygon = lower_polygon
            for higher_index, (higher_name, higher_polygon) in reversed(list(enumerate(polygon_dict.items()))[:lower_index]):
                diff_polygon = diff_polygon.difference(higher_polygon)
            polygons_tiled_dict[lower_name] = diff_polygon

        # print("Tiled")
        # for key in polygons_tiled_dict.keys():
        #     print(key, polygons_tiled_dict[key])

        # Break up polygon edges so that plane is tiled with no partially overlapping line segments
        polygons_broken_dict = {}
        for first_index, (first_name, first_polygons) in enumerate(polygon_dict.items()):
            first_polygons = polygons_tiled_dict[first_name]
            if first_polygons.type == "MultiPolygon":
                multi = True
                first_polygons = first_polygons.geoms
            else:
                multi= False
                first_polygons = [first_polygons]
            broken_polygons = []
            for first_polygon in first_polygons:
                # Exterior
                first_exterior_line = LineString(first_polygon.exterior)
                for second_index, (second_name, second_polygons) in enumerate(polygon_dict.items()):
                    if second_name == first_name:
                        continue
                    else:
                        second_polygons = polygons_tiled_dict[first_name]
                        if second_polygons.type == "MultiPolygon":
                            second_polygons = second_polygons.geoms
                        else:
                            second_polygons = [second_polygons]
                        for second_polygon in second_polygons:
                            # Exterior
                            second_exterior_line = LineString(second_polygon.exterior)
                            intersections = first_exterior_line.intersection(second_exterior_line)
                            if intersections.is_empty:
                                continue
                            else:
                                for intersection in intersections.geoms:
                                    new_coords_start, new_coords_end = intersection.boundary.geoms
                                    first_exterior_line = linemerge(split(first_exterior_line, new_coords_start))
                                    first_exterior_line = linemerge(split(first_exterior_line, new_coords_end))
                            # Interiors
                            for second_interior_line in second_polygon.interiors:
                                second_interior_line = LineString(second_interior_line)
                                intersections = first_exterior_line.intersection(second_interior_line)
                                if intersections.is_empty:
                                    continue
                                else:
                                    for intersection in intersections.geoms:
                                        new_coords_start, new_coords_end = intersection.boundary.geoms
                                        first_exterior_line = linemerge(split(first_exterior_line, new_coords_start))
                                        first_exterior_line = linemerge(split(first_exterior_line, new_coords_end))
                # Interiors
                first_polygon_interiors = []
                for first_interior_line in first_polygon.interiors:
                    for second_index, (second_name, second_polygons) in enumerate(polygon_dict.items()):
                        if second_name == first_name:
                            continue
                        else:
                            second_polygons = polygons_tiled_dict[first_name]
                            if second_polygons.type == "MultiPolygon":
                                second_polygons = second_polygons.geoms
                            else:
                                second_polygons = [second_polygons]
                            for second_polygon in second_polygons:
                                # Exterior
                                second_exterior_line = LineString(second_polygon.exterior)
                                intersections = first_interior_line.intersection(second_exterior_line)
                                if intersections.is_empty:
                                    continue
                                else:
                                    for intersection in intersections.geoms:
                                        new_coords_start, new_coords_end = intersection.boundary.geoms
                                        first_interior_line = linemerge(split(first_interior_line, new_coords_start))
                                        first_interior_line = linemerge(split(first_interior_line, new_coords_end))
                                # Interiors
                                for second_interior_line in second_polygon.interiors:
                                    second_interior_line = LineString(second_interior_line)
                                    intersections = first_interior_line.intersection(second_interior_line)
                                    if intersections.is_empty:
                                        continue
                                    else:
                                        for intersection in intersections.geoms:
                                            new_coords_start, new_coords_end = intersection.boundary.geoms
                                            first_interior_line = linemerge(split(first_interior_line, new_coords_start))
                                            first_interior_line = linemerge(split(first_interior_line, new_coords_end))
                    first_polygon_interiors.append(first_interior_line)
                broken_polygons.append(Polygon(first_exterior_line, holes=first_polygon_interiors))
            if multi == True:
                polygons_broken_dict[first_name] = MultiPolygon(broken_polygons)
            else:
                polygons_broken_dict[first_name] = broken_polygons[0]

        print("Broken")
        for key in polygons_broken_dict.keys():
            print(key, polygons_broken_dict[key])
        
        # Add surfaces, reusing lines to simplify at early stage
        meshtracker = MeshTracker(model=model)
        for polygon_name, polygon in reversed(polygons_broken_dict.items()):
            plane_surface = meshtracker.add_xy_surface(polygon, polygon_name)
            model.add_physical(plane_surface, f"{polygon_name}")

        # Refinement in surfaces
        n = 0
        refinement_fields = []
        for label, resolution in resolutions.items():
            # Inside surface
            mesh_resolution = resolution["resolution"]
            gmsh.model.mesh.field.add("MathEval", n)
            gmsh.model.mesh.field.setString(n, "F", f"{mesh_resolution}")
            gmsh.model.mesh.field.add("Restrict", n+1)
            gmsh.model.mesh.field.setNumber(n+1, "InField", n)
            gmsh.model.mesh.field.setNumbers(n+1, "SurfacesList", meshtracker.get_gmsh_xy_surfaces_from_label(label))
            # Around surface
            mesh_distance = resolution["distance"]
            gmsh.model.mesh.field.add("Distance", n+2)
            gmsh.model.mesh.field.setNumbers(n+2, "CurvesList", meshtracker.get_gmsh_xy_lines_from_label(label))
            gmsh.model.mesh.field.setNumber(n+2, "Sampling", 100)
            gmsh.model.mesh.field.add("Threshold", n+3)
            gmsh.model.mesh.field.setNumber(n+3, "InField", n+2)
            gmsh.model.mesh.field.setNumber(n+3, "SizeMin", mesh_resolution)
            gmsh.model.mesh.field.setNumber(n+3, "SizeMax", default_resolution_max)
            gmsh.model.mesh.field.setNumber(n+3, "DistMin", 0)
            gmsh.model.mesh.field.setNumber(n+3, "DistMax", mesh_distance)
            # Save and increment
            refinement_fields.append(n+1)
            refinement_fields.append(n+3)
            n += 4

        # Use the smallest element size overall
        gmsh.model.mesh.field.add("Min", n)
        gmsh.model.mesh.field.setNumbers(n, "FieldsList", refinement_fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(n)

        gmsh.model.mesh.MeshSizeFromPoints = 0
        gmsh.model.mesh.MeshSizeFromCurvature = 0
        gmsh.model.mesh.MeshSizeExtendFromBoundary = 0

        # Fuse edges (bandaid)
        # gmsh.model.occ.synchronize()
        # gmsh.model.occ.removeAllDuplicates()
        # gmsh.model.occ.synchronize()

        # Extract all unique lines (TODO: identify interfaces in label)
        i = 0
        for line in meshtracker.gmsh_xy_lines:
            model.add_physical(line, f"line_{i}")
            i += 1

        if filename:
            mesh = geometry.generate_mesh(dim=2, verbose=True)
            gmsh.write(f"{filename}")
        else:
            mesn = geometry.generate_mesh(dim=2, verbose=True)

        return mesh
        

if __name__ == "__main__":

    import gmsh

    wsim = 2
    hclad = 2
    hbox = 2
    offset_core = -0.1
    offset_core2 = 1
    wcore = 0.5
    hcore = 0.22
    core = Polygon([
            Point(-wcore/2, -hcore/2 + offset_core),
            Point(-wcore/2, hcore/2 + offset_core),
            Point(wcore/2, hcore/2 + offset_core),
            Point(wcore/2, -hcore/2 + offset_core),
        ])
    core2 = Polygon([
            Point(-wcore/2, -hcore/2 + offset_core2),
            Point(-wcore/2, hcore/2 + offset_core2),
            Point(wcore/2, hcore/2 + offset_core2),
            Point(wcore/2, -hcore/2 + offset_core2),
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

    polygons = OrderedDict()
    polygons["core"] = core 
    polygons["core2"] = core2
    polygons["clad"] = clad
    polygons["box"] = box

    resolutions = {}
    resolutions["core"] = {"resolution": 0.01, "distance": 5}
    resolutions["core2"] = {"resolution": 0.01, "distance": 5}
    # resolutions["clad"] = {"resolution": 0.1, "dist_min": 0.01, "dist_max": 0.3}


    mesh = mesh_from_polygons(polygons, resolutions, filename="mesh.msh")

    # gmsh.write("mesh.msh")
    # gmsh.clear()
    # mesh.__exit__()

    import meshio

    mesh_from_file = meshio.read("mesh.msh")

    def create_mesh(mesh, cell_type, prune_z=True):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        points = mesh.points
        return meshio.Mesh(
            points=points,
            cells={cell_type: cells},
            cell_data={"name_to_read": [cell_data]},
        )

    line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
    meshio.write("facet_mesh.xdmf", line_mesh)

    triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
    meshio.write("mesh.xdmf", triangle_mesh)