from collections import OrderedDict
from itertools import combinations, product
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
from shapely.ops import linemerge, polygonize, split, unary_union

from femwell.mesh.meshtracker import MeshTracker


def break_line_(line, other_line):
    intersections = line.intersection(other_line)
    if not intersections.is_empty:
        for intersection in (
            intersections.geoms if hasattr(intersections, "geoms") else [intersections]
        ):
            # if type == "", intersection.type != 'Point':
            if intersection.geom_type == "Point":
                line = linemerge(split(line, intersection))
            else:
                new_coords_start, new_coords_end = intersection.boundary.geoms
                line = linemerge(split(line, new_coords_start))
                line = linemerge(split(line, new_coords_end))
    return line


def mesh_from_Dict(
    shapes_dict: Dict,
    resolutions: Optional[Dict[str, Dict[str, float]]] = None,
    default_resolution_min: float = 0.01,
    default_resolution_max: float = 0.5,
    filename: Optional[str] = None,
    gmsh_algorithm: int = 5,
    global_quad: Optional[bool] = False,
    verbose: bool = False,
):
    """
    Given a dict of shapely Polygons, creates a mesh conformal to all the BooleanFragment surfaces and interfaces from the intersection of all polygons.
    Returns a mesh with physicals corresponding to the original shapes_dict boundaries (which may comprise multiple surface entities)
    """

    with pygmsh.occ.geometry.Geometry() as geometry:
        # geometry = pygmsh.occ.geometry.Geometry()
        geometry.characteristic_length_min = default_resolution_min
        geometry.characteristic_length_max = default_resolution_max
        gmsh.option.setNumber("Mesh.Algorithm", gmsh_algorithm)

        model = geometry

        # Add overall bounding box
        total = unary_union(list(shapes_dict.values()))
        all_polygons = [total, *list(shapes_dict.values())]
        # Break up all shapes so that plane is tiled with non-overlapping layers, get the maximal number of fragments
        # Equivalent to BooleanFragments
        listpoly = [a.intersection(b) for a, b in combinations(all_polygons, 2)]
        rings = [
            LineString(list(object.exterior.coords))
            for object in listpoly
            if not (object.geom_type in ["Point", "LineString"] or object.is_empty)
        ]

        union = unary_union(rings)
        shapes_tiled_list = list(polygonize(union))

        # Add surfaces, logging lines
        meshtracker = MeshTracker(model=model)
        for i, polygon in enumerate(shapes_tiled_list):
            meshtracker.add_xy_surface(polygon, i, physical=False)
        # Tag physicals
        # TODO integrate in meshtracker
        surface_name_mapping = {}
        surface_resolution_mapping = {}
        for polygon_name, polygons in shapes_dict.items():
            surfaces = []
            surface_indices = []
            for polygon in polygons.geoms if hasattr(polygons, "geoms") else [polygons]:
                # Identify enclosed surfaces
                for index in range(len(meshtracker.shapely_xy_surfaces)):
                    if polygon.contains(meshtracker.shapely_xy_surfaces[index]):
                        surfaces.append(meshtracker.gmsh_xy_surfaces[index])
                        surface_indices.append(index)
                        if resolutions:
                            if index in surface_resolution_mapping:
                                surface_resolution_mapping[index] = min(
                                    resolutions[polygon_name]["resolution"],
                                    surface_resolution_mapping[index],
                                )
                            else:
                                surface_resolution_mapping[index] = min(
                                    resolutions[polygon_name]["resolution"],
                                    default_resolution_max,
                                )
            meshtracker.model.add_physical(surfaces, polygon_name)
            surface_name_mapping[polygon_name] = surface_indices

        current_resolutions = []
        if resolutions:
            # Refinement in surfaces
            n = 0
            refinement_fields = []
            for surface_indices in surface_name_mapping.values():
                for surface_index in surface_indices:
                    gmsh.model.mesh.field.add("MathEval", n)
                    gmsh.model.mesh.field.setString(
                        n, "F", f"{surface_resolution_mapping[surface_index]}"
                    )
                    gmsh.model.mesh.field.add("Restrict", n + 1)
                    gmsh.model.mesh.field.setNumber(n + 1, "InField", n)
                    gmsh.model.mesh.field.setNumbers(
                        n + 1,
                        "SurfacesList",
                        [meshtracker.gmsh_xy_surfaces[surface_index]._id],
                    )
                    # # Around surface
                    # mesh_distance = mesh_setting["distance"]
                    # gmsh.model.mesh.field.add("Distance", n+2)
                    # gmsh.model.mesh.field.setNumbers(n+2, "CurvesList", meshtracker.get_gmsh_xy_lines_from_label(label))
                    # gmsh.model.mesh.field.setNumber(n+2, "Sampling", 100)
                    # gmsh.model.mesh.field.add("Threshold", n+3)
                    # gmsh.model.mesh.field.setNumber(n+3, "InField", n+2)
                    # gmsh.model.mesh.field.setNumber(n+3, "SizeMin", mesh_resolution)
                    # gmsh.model.mesh.field.setNumber(n+3, "SizeMax", default_resolution_max)
                    # gmsh.model.mesh.field.setNumber(n+3, "DistMin", 0)
                    # gmsh.model.mesh.field.setNumber(n+3, "DistMax", mesh_distance)
                    # Save and increment
                    refinement_fields.append(n + 1)
                    # refinement_fields.append(n+3)
                    n += 2

            if global_quad:
                gmsh.option.setNumber("Mesh.Algorithm", 8)
                gmsh.option.setNumber("Mesh.RecombineAll", 1)

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
        # i = 0
        # for index, line in enumerate(meshtracker.gmsh_xy_segments):
        #     model.add_physical(line, f"{meshtracker.xy_segments_main_labels[index]}_{meshtracker.xy_segments_secondary_labels[index]}_{i}")
        #     i += 1

        # For periodicity

        mesh = geometry.generate_mesh(dim=2, verbose=verbose)

        if filename:
            gmsh.write(f"{filename}")

        return mesh


def mesh_from_OrderedDict(
    shapes_dict: OrderedDict,
    resolutions: Optional[Dict[str, Dict[str, float]]] = None,
    default_resolution_min: float = 0.01,
    default_resolution_max: float = 0.5,
    filename: Optional[str] = None,
    gmsh_algorithm: int = 5,
    global_quad: Optional[bool] = False,
    verbose: bool = False,
    periodic_lines: Optional[Tuple[(str, str)]] = None,
    mesh_scaling_factor: float = 1.0,
):
    """
    Given an ordered dict of shapely Polygons, creates a mesh containing polygon surfaces according to the dict order.
    Returns a gmsh msh with physicals corresponding to the shapes_dict boundaries (which is the minimal number of surfaces for each key)

    periodic_lines: (label1, label1) tuples forcing the mesh of line[label1] to map to the mesh of line[label2]. Currently only works if the lines are not intersected.
    """

    with pygmsh.occ.geometry.Geometry() as geometry:
        # geometry = pygmsh.occ.geometry.Geometry()
        geometry.characteristic_length_min = default_resolution_min
        geometry.characteristic_length_max = default_resolution_max
        gmsh.option.setNumber("Mesh.Algorithm", gmsh_algorithm)

        model = geometry

        # Break up shapes in order so that plane is tiled with non-overlapping layers, overriding shapes according to an order
        shapes_tiled_dict = OrderedDict()
        for lower_index, (lower_name, lower_shape) in reversed(
            list(enumerate(shapes_dict.items()))
        ):
            diff_shape = lower_shape
            for higher_index, (higher_name, higher_shape) in reversed(
                list(enumerate(shapes_dict.items()))[:lower_index]
            ):
                diff_shape = diff_shape.difference(higher_shape)
            shapes_tiled_dict[lower_name] = diff_shape

        # Break up lines and polygon edges so that plane is tiled with no partially overlapping line segments
        polygons_broken_dict = OrderedDict()
        lines_broken_dict = OrderedDict()
        for first_index, (first_name, first_shape) in enumerate(shapes_dict.items()):
            first_shape = shapes_tiled_dict[first_name]
            broken_shapes = []
            for first_shape in (
                first_shape.geoms if hasattr(first_shape, "geoms") else [first_shape]
            ):
                # First line exterior
                first_exterior_line = (
                    LineString(first_shape.exterior)
                    if first_shape.geom_type == "Polygon"
                    else first_shape
                )
                for second_index, (second_name, second_shapes) in enumerate(shapes_dict.items()):
                    # Do not compare to itself
                    if second_name == first_name:
                        continue
                    else:
                        second_shapes = shapes_tiled_dict[second_name]
                        for second_shape in (
                            second_shapes.geoms
                            if hasattr(second_shapes, "geoms")
                            else [second_shapes]
                        ):
                            # Second line exterior
                            second_exterior_line = (
                                LineString(second_shape.exterior)
                                if second_shape.geom_type == "Polygon"
                                else second_shape
                            )
                            first_exterior_line = break_line_(
                                first_exterior_line, second_exterior_line
                            )
                            # Second line interiors
                            for second_interior_line in (
                                second_shape.interiors
                                if second_shape.geom_type == "Polygon"
                                else []
                            ):
                                second_interior_line = LineString(second_interior_line)
                                first_exterior_line = break_line_(
                                    first_exterior_line, second_interior_line
                                )
                # First line interiors
                if first_shape.geom_type in ["Polygon", "MultiPolygon"]:
                    first_shape_interiors = []
                    for first_interior_line in first_shape.interiors:
                        first_interior_line = LineString(first_interior_line)
                        for second_index, (second_name, second_shapes) in enumerate(
                            shapes_dict.items()
                        ):
                            if second_name == first_name:
                                continue
                            else:
                                second_shapes = shapes_tiled_dict[second_name]
                                for second_shape in (
                                    second_shapes.geoms
                                    if hasattr(second_shapes, "geoms")
                                    else [second_shapes]
                                ):
                                    # Exterior
                                    second_exterior_line = (
                                        LineString(second_shape.exterior)
                                        if second_shape.geom_type == "Polygon"
                                        else second_shape
                                    )
                                    first_interior_line = break_line_(
                                        first_interior_line, second_exterior_line
                                    )
                                    # Interiors
                                    for second_interior_line in (
                                        second_shape.interiors
                                        if second_shape.geom_type == "Polygon"
                                        else []
                                    ):
                                        second_interior_line = LineString(second_interior_line)
                                        intersections = first_interior_line.intersection(
                                            second_interior_line
                                        )
                                        first_interior_line = break_line_(
                                            first_interior_line, second_interior_line
                                        )
                        first_shape_interiors.append(first_interior_line)
                if first_shape.geom_type in ["Polygon", "MultiPolygon"]:
                    broken_shapes.append(Polygon(first_exterior_line, holes=first_shape_interiors))
                else:
                    broken_shapes.append(LineString(first_exterior_line))
            if first_shape.geom_type in ["Polygon", "MultiPolygon"]:
                polygons_broken_dict[first_name] = (
                    MultiPolygon(broken_shapes) if len(broken_shapes) > 1 else broken_shapes[0]
                )
            else:
                lines_broken_dict[first_name] = (
                    MultiLineString(broken_shapes) if len(broken_shapes) > 1 else broken_shapes[0]
                )

        # Add lines, reusing line segments
        meshtracker = MeshTracker(model=model)
        for line_name, line in lines_broken_dict.items():
            meshtracker.add_get_xy_line(line, line_name)

        # Add surfaces, reusing lines
        for polygon_name, polygon in polygons_broken_dict.items():
            meshtracker.add_xy_surface(polygon, polygon_name)

        # Embed lines in surfaces if required
        for index_surface in range(len(meshtracker.shapely_xy_surfaces)):
            polygon = meshtracker.shapely_xy_surfaces[index_surface]
            for index_segment in range(len(meshtracker.shapely_xy_segments)):
                line = meshtracker.shapely_xy_segments[index_segment]

                intersection = line - polygon.exterior
                if not intersection.is_empty and polygon.contains(intersection):
                    model.in_surface(
                        meshtracker.gmsh_xy_segments[index_segment],
                        meshtracker.gmsh_xy_surfaces[index_surface],
                    )

        # Refinement in surfaces
        n = 0
        refinement_fields = []
        for label, mesh_setting in resolutions.items():
            # Inside surface
            mesh_resolution = mesh_setting["resolution"]
            gmsh.model.mesh.field.add("MathEval", n)
            gmsh.model.mesh.field.setString(n, "F", f"{mesh_resolution}")
            gmsh.model.mesh.field.add("Restrict", n + 1)
            gmsh.model.mesh.field.setNumber(n + 1, "InField", n)
            gmsh.model.mesh.field.setNumbers(
                n + 1,
                "SurfacesList",
                meshtracker.get_gmsh_xy_surfaces_from_label(label),
            )
            # Around surface
            mesh_distance = mesh_setting["distance"]
            gmsh.model.mesh.field.add("Distance", n + 2)
            gmsh.model.mesh.field.setNumbers(
                n + 2, "CurvesList", meshtracker.get_gmsh_xy_lines_from_label(label)
            )
            gmsh.model.mesh.field.setNumber(n + 2, "Sampling", 100)
            gmsh.model.mesh.field.add("Threshold", n + 3)
            gmsh.model.mesh.field.setNumber(n + 3, "InField", n + 2)
            gmsh.model.mesh.field.setNumber(n + 3, "SizeMin", mesh_resolution)
            gmsh.model.mesh.field.setNumber(n + 3, "SizeMax", default_resolution_max)
            gmsh.model.mesh.field.setNumber(n + 3, "DistMin", 0)
            gmsh.model.mesh.field.setNumber(n + 3, "DistMax", mesh_distance)
            # Save and increment
            refinement_fields.append(n + 1)
            refinement_fields.append(n + 3)
            n += 4

        if global_quad:
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.option.setNumber("Mesh.RecombineAll", 1)

        # Use the smallest element size overall
        gmsh.model.mesh.field.add("Min", n)
        gmsh.model.mesh.field.setNumbers(n, "FieldsList", refinement_fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(n)

        gmsh.model.mesh.MeshSizeFromPoints = 0
        gmsh.model.mesh.MeshSizeFromCurvature = 0
        gmsh.model.mesh.MeshSizeExtendFromBoundary = 0

        # Tag all interfacial lines
        for surface1, surface2 in combinations(polygons_broken_dict.keys(), 2):
            interfaces = []
            for index, line in enumerate(meshtracker.gmsh_xy_segments):
                if (
                    meshtracker.xy_segments_main_labels[index] == surface1
                    and meshtracker.xy_segments_secondary_labels[index] == surface2
                ) or (
                    meshtracker.xy_segments_main_labels[index] == surface2
                    and meshtracker.xy_segments_secondary_labels[index] == surface1
                ):
                    interfaces.append(line)
            if interfaces:
                model.add_physical(interfaces, f"{surface1}___{surface2}")

        gmsh.model.occ.synchronize()

        # Force periodicity (experimental)
        def validate_lines(line1, line2):
            """TODO create a module for validating geometries."""

        if periodic_lines:
            for label1, label2 in periodic_lines:
                # if validate_lines(): # TODO
                line1 = shapes_dict[label1]
                line2 = shapes_dict[label2]
                gmsh.model.setCurrent("pygmsh model")
                translation = np.array(line1.coords[0]) - np.array(line2.coords[0])
                gmsh.model.mesh.setPeriodic(
                    1,
                    meshtracker.get_gmsh_xy_lines_from_label(label1),
                    meshtracker.get_gmsh_xy_lines_from_label(label2),
                    [1, 0, 0, translation[0], 0, 1, 0, translation[1], 0, 0, 1, 0, 0, 0, 0, 1],
                )
                # else: # TODO
                #     raise ValueError("Periodic line pairs must be parallel and have the same straight length in the final, intersected geometry.")

        gmsh.option.setNumber("Mesh.ScalingFactor", mesh_scaling_factor)
        mesh = geometry.generate_mesh(dim=2, verbose=verbose)

        if filename:
            gmsh.write(f"{filename}")

        import contextlib
        import tempfile

        import meshio

        with contextlib.redirect_stdout(None):
            with tempfile.TemporaryDirectory() as tmpdirname:
                gmsh.write(f"{tmpdirname}/mesh.msh")
                return meshio.read(f"{tmpdirname}/mesh.msh")


if __name__ == "__main__":
    from collections import OrderedDict

    import gmsh
    from mesh import mesh_from_OrderedDict
    from shapely.geometry import LineString, Polygon

    width = 4
    length = 10.5
    pml = 0.5

    width_wg_1 = 0.5
    length_wg_1 = 5

    width_wg_2 = 2
    length_wg_2 = 5

    core = Polygon(
        [
            (-width_wg_1 / 2, -length_wg_1),
            (-width_wg_1 / 2, 0),
            (-width_wg_2 / 2, 0),
            (-width_wg_2 / 2, length_wg_2),
            (width_wg_2 / 2, length_wg_2),
            (width_wg_2 / 2, 0),
            (width_wg_1 / 2, 0),
            (width_wg_1 / 2, -length_wg_1),
        ]
    )
    print(core)

    box = Polygon(
        [
            (-width_wg_2 - 1, -6),
            (-width_wg_2 - 1, 6),
            (width_wg_2 - 0.5, 6),
            (width_wg_2 - 0.5, -6),
        ]
    )
    print(box)

    source = LineString([(width_wg_2 / 2, -length_wg_1 / 2), (-width_wg_2 / 2, -length_wg_1 / 2)])
    print(source)

    left_wall_up = LineString([(-width_wg_2 - 1, -2), (-width_wg_2 - 1, 6)])
    right_wall_up = LineString([(width_wg_2 - 0.5, -2), (width_wg_2 - 0.5, 6)])
    left_wall_dw = LineString([(-width_wg_2 - 1, -6), (-width_wg_2 - 1, -2)])
    right_wall_dw = LineString([(width_wg_2 - 0.5, -6), (width_wg_2 - 0.5, -2)])

    polygons = OrderedDict(
        left_wall_up=left_wall_up,
        right_wall_up=right_wall_up,
        left_wall_dw=left_wall_dw,
        right_wall_dw=right_wall_dw,
        source=source,
        core=core,
        box=box,
        # pml=core.buffer(2, resolution=4) - core.buffer(1, resolution=4),
    )

    resolutions = dict(
        source={"resolution": 0.02, "distance": 1},
        core={"resolution": 0.02, "distance": 1},
        # left_wall_up={"resolution": 1, "distance": 1},
        # right_wall_up={"resolution": 0.5, "distance": 1},
    )

    mesh = mesh_from_OrderedDict(
        polygons,
        resolutions,
        filename="mesh.msh",
        default_resolution_max=1,
        # periodic_lines=[("left_wall_up", "right_wall_up"), ("left_wall_dw", "right_wall_dw")],
        mesh_scaling_factor=1e-4,
    )

    # wmode = 1
    # wsim = 2
    # hclad = 2
    # hbox = 2
    # wcore = 0.5
    # hcore = 0.22
    # offset_core = -0.1
    # offset_core2 = 1

    # # Lines can be added, which is useful to define boundary conditions at various simulation edges
    # left_edge = LineString([Point(-wsim/2, -hcore/2  - hbox),
    #                         Point(-wsim/2, -hcore/2 + hclad)])
    # right_edge = LineString([Point(wsim/2, -hcore/2  - hbox),
    #                         Point(wsim/2, -hcore/2 + hclad)])
    # top_edge = LineString([Point(-wsim/2, -hcore/2 + hclad),
    #                         Point(wsim/2, -hcore/2 + hclad)])
    # bottom_edge = LineString([Point(-wsim/2, -hcore/2  - hbox),
    #                         Point(wsim/2, -hcore/2  - hbox)])

    # # Polygons not only have an edge, but an interior
    # core = Polygon([
    #         Point(-wcore/2, -hcore/2 + offset_core),
    #         Point(-wcore/2, hcore/2 + offset_core),
    #         Point(wcore/2, hcore/2 + offset_core),
    #         Point(wcore/2, -hcore/2 + offset_core),
    #     ])
    # core2 = Polygon([
    #         Point(-wcore/2, -hcore/2 + offset_core2),
    #         Point(-wcore/2, hcore/2 + offset_core2),
    #         Point(wcore/2, hcore/2 + offset_core2),
    #         Point(wcore/2, -hcore/2 + offset_core2),
    #     ])
    # clad = Polygon([
    #         Point(-wsim/2, -hcore/2),
    #         Point(-wsim/2, -hcore/2 + hclad),
    #         Point(wsim/2, -hcore/2 + hclad),
    #         Point(wsim/2, -hcore/2),
    #     ])
    # box = Polygon([
    #         Point(-wsim/2, -hcore/2),
    #         Point(-wsim/2, -hcore/2 - hbox),
    #         Point(wsim/2, -hcore/2 - hbox),
    #         Point(wsim/2, -hcore/2),
    #     ])

    # # The order in which objects are inserted into the OrderedDict determines overrrides
    # # shapes = OrderedDict()
    # shapes = {}
    # # shapes["left_edge"] = left_edge
    # # shapes["right_edge"] = right_edge
    # # shapes["top_edge"] = top_edge
    # # shapes["bottom_edge"] = bottom_edge
    # shapes["core"] = core
    # shapes["core2"] = core2
    # shapes["clad"] = clad
    # shapes["box"] = box

    # # The resolution dict is not ordered, and can be used to set mesh resolution at various element
    # # The edge of a polygon and another polygon (or entire simulation domain) will form a line object that can be refined independently
    # resolutions = {}
    # resolutions["core"] = {"resolution": 0.02, "distance": 2}
    # resolutions["core2"] = {"resolution": 0.02, "distance": 2}
    # resolutions["clad"] = {"resolution": 0.5, "distance": 2}
    # resolutions["box"] = {"resolution": 0.5, "distance": 2}
    # # resolutions["core_clad"] = {"resolution": 0.05, "distance": 0.5}
    # # resolutions["clad_box"] = {"resolution": 0.05, "distance": 0.5}
    # # resolutions["bottom_edge"] = {"resolution": 0.05, "distance": 0.5}
    # # resolutions["left_edge"] = {"resolution": 0.05, "distance": 0.5}
    # # resolutions["clad"] = {"resolution": 0.1, "dist_min": 0.01, "dist_max": 0.3}

    # quad = False
    # # mesh = mesh_from_OrderedDict(shapes, resolutions, filename="mesh.msh", default_resolution_max=.3, global_quad=quad)
    # mesh = mesh_from_Dict(shapes, resolutions, filename="mesh.msh", default_resolution_max=0.5, global_quad=quad)

    # # gmsh.write("mesh.msh")
    # # gmsh.clear()
    # # mesh.__exit__()

    # import meshio

    # mesh_from_file = meshio.read("mesh.msh")

    # def create_mesh(mesh, cell_type, prune_z=True):
    #     cells = mesh.get_cells_type(cell_type)
    #     cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    #     points = mesh.points
    #     return meshio.Mesh(
    #         points=points,
    #         cells={cell_type: cells},
    #         cell_data={"name_to_read": [cell_data]},
    #     )

    # # line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
    # # meshio.write("facet_mesh.xdmf", line_mesh)

    # if quad == True:
    #     triangle_mesh = create_mesh(mesh_from_file, "quad", prune_z=True)
    # else:
    #     triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
    # meshio.write("mesh.xdmf", triangle_mesh)
