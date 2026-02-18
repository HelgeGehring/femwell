from collections import OrderedDict
from typing import Dict, Optional

import gdsfactory as gf
import numpy as np
from gdsfactory.component import Component
from gdsfactory.simulation.gmsh import (
    cleanup_component,
    get_uz_bounds_layers,
    order_layerstack,
)
from gdsfactory.tech import LayerStack
from shapely.affinity import translate
from shapely.geometry import LineString, MultiPoint, MultiPolygon, Polygon

# from femwell.mesh import mesh_from_Dict

nm = 1e-3


def to_polygons(geometries):
    for geometry in geometries:
        if isinstance(geometry, Polygon):
            yield geometry
        else:
            yield from geometry


def get_vertices(polygon):
    """Return all polygon vertices (interior and exterior)"""
    vertices = []
    for polygon_hole in list(polygon.interiors):
        vertices.extend(iter(MultiPoint(polygon_hole.coords).geoms))
    # Parse boundary
    vertices.extend(iter(MultiPoint(polygon.exterior.coords).geoms))

    return vertices


def get_polygon_x_bounds(polygon):
    """Return x_bounds of polygon vertices.

    Propagation direction is "x" in component ("z" in mode solver)
    """
    return [p.x for p in get_vertices(polygon)]


def get_unique_x_bounds(layer_polygons_dict):
    """Return unique x_bounds across all polygon vertices of a layer_polygons_dict

    Propagation direction is "x" in component ("z" in mode solver)
    """
    xs = [0.0]
    for polygons in layer_polygons_dict.values():
        for polygon in polygons.geoms if hasattr(polygons, "geoms") else [polygons]:
            xs.extend(get_polygon_x_bounds((polygon)))
    return np.sort(np.unique(np.array(xs)))


def get_mode_regions(component, layerstack, tol=1e-6):
    """Return x_bounds of polygon vertices where there is a change of cross-section.

    Propagation direction is "x" in component ("z" in mode solver)

    Args:
        component: gdsfactory Component to process
        layerstack: gdsfactory LayerStack to process
    Returns:
        x_changing_bounds: x-bounds where the component needs to be meshed
        x_not_changing_bounds: x-bounds where there the structure is unchanging (free propagation)
    """
    bbox = gf.components.bbox(bbox=component.bbox).get_polygons()[0][:, 1]
    ymin = np.min(bbox)
    ymax = np.max(bbox)
    layer_polygons_dict = cleanup_component(component, layerstack)
    x_bounds = get_unique_x_bounds(layer_polygons_dict)

    x_changing_bounds = []
    x_not_changing_bounds = []
    for x1, x2 in [(x_bounds[i] + tol, x_bounds[i + 1] - tol) for i in range(len(x_bounds) - 1)]:
        found_different = False
        for layername, polygons in layer_polygons_dict.items():
            if found_different:
                continue
            line_x1 = LineString([[x1, ymin], [x1, ymax]])
            line_x2 = LineString([[x2, ymin], [x2, ymax]])

            xsection_x1 = polygons.intersection(line_x1)
            xsection_x2 = polygons.intersection(line_x2)

            if not xsection_x1.equals(translate(xsection_x2, xoff=x1 - x2)):
                found_different = True
                x_changing_bounds.append([x1 - tol, x2 + tol])

        if not found_different:
            x_not_changing_bounds.append([x1 - tol, x2 + tol])

    return x_changing_bounds, x_not_changing_bounds


def slice_component_xbounds(component, layerstack, mesh_step=100 * nm):
    """Returns minimal list of x-coordinates where cross-section is to be taken."""
    # Process polygon to extract regions of change and free propagation
    x_changing_bounds, x_not_changing_bounds = get_mode_regions(component, layerstack, tol=1e-6)

    # Where geometry is changing, mesh according to mesh_step
    x_coordinates = [np.arange(x1, x2, mesh_step) for x1, x2 in x_changing_bounds]
    # Where not changing, just return the bounds
    x_coordinates.extend(np.array([x1]) for x1, x2 in x_not_changing_bounds)
    # Return sorted bounds
    return np.sort(np.concatenate(x_coordinates).ravel())


"""Below fails due to complicated mesh to generate."""

# def slice_component_polygons(component, layerstack, ymin=-10, ymax=10, mesh_step=100*nm):
#     """Returns a dict of dicts "layer__x" of polygons.

#         "layer__x" to be used as the gmsh physical label.
#     """
#     layer_polygons_dict = cleanup_component(component, layerstack)
#     x_coords = slice_component_xbounds(component, layerstack, mesh_step=100*nm)
#     layer_order = order_layerstack(layerstack)
#     shapes = {}
#     for x_coord in x_coords:
#         xsection_bounds = [[x_coord, ymin],[x_coord, ymax]]
#         bounds_dict = get_uz_bounds_layers(layer_polygons_dict, xsection_bounds, layerstack)
#         for layer in layer_order:
#             layer_shapes = []
#             for polygon in bounds_dict[layer]:
#                 layer_shapes.append(polygon)
#             shapes[f"{layer}_{x_coord}"] = MultiPolygon(to_polygons(layer_shapes))
#     return shapes


# def mesh_from_slices(component: Component,
#                         layerstack: LayerStack,
#                         mesh_step: float = 500 * nm,
#                         resolutions: Optional[Dict[str, Dict[str, float]]] = {},
#                         default_resolution_min: float = 0.01,
#                         default_resolution_max: float = 0.5,
#                         filename: Optional[str] = None,
#                         gmsh_algorithm: int = 5,
#                         global_quad: Optional[bool] = False,
#                     ):
#     """Returns a single cross-sectional uz mesh conditioned on multiple slices."""
#     shapes = slice_component_polygons(component, layerstack, mesh_step)
#     return mesh_from_Dict(
#             shapes_dict=shapes,
#             resolutions=resolutions,
#             filename='test.msh'
#         )


if __name__ == "__main__":
    import gdsfactory as gf
    import gdsfactory.simulation.gmsh as gfmesh

    """
    The below mininmal example shows how to
    (1) Define a gsfactory component
    (2) Define a reduced LayerStack
    (3) Find the x-coordinates where there are interfaces
    (4) Acquire a mesh at a given x-coordinate of the component
    """

    """
    (1) Get component
    """
    c = gf.Component()

    taper = c.add_ref(
        gf.components.taper(
            length=10.0,
            width1=0.5,
            width2=2,
            cross_section="rib",
        )
    )
    straight = c.add_ref(gf.components.straight(10, width=2, cross_section="rib"))
    straight.connect("o1", taper.ports["o2"])
    straight2 = c.add_ref(gf.components.straight(10, width=2, cross_section="strip"))
    straight2.connect("o1", straight.ports["o2"])
    straight = c.add_ref(gf.components.straight(10, width=2, cross_section="rib"))
    straight.connect("o1", straight2.ports["o2"])
    taper = c.add_ref(
        gf.components.taper(
            length=10.0,
            width1=0.5,
            width2=2,
            cross_section="rib",
        )
    )
    taper.connect("o2", straight.ports["o2"])
    bend = c.add_ref(
        gf.components.bend_euler(
            angle=20.0,
            p=0.5,
            width=0.5,
            cross_section="strip",
        ).move([12, -5])
    )

    """
    (2) Get LayerStack
    """
    import numpy as np
    from gdsfactory.tech import LayerStack, get_layer_stack_generic

    filtered_layerstack = LayerStack(
        layers={
            k: get_layer_stack_generic().layers[k]
            for k in (
                "slab90",
                "core",
            )
        }
    )

    """
    (3) Get x-coordinates
    """
    # Compute the x-coordinates where the cross-section changes
    lines_x = slice_component_xbounds(c, filtered_layerstack, mesh_step=500 * nm)
    for x in lines_x:
        P = gf.Path([[x, -20], [x, 20]])
        X = gf.CrossSection(width=0.001, layer=(99, 0))
        line = gf.path.extrude(P, X)
        c << line

    c.show()

    """
    (4) Get a simulatoin a x
    """
    resolutions = {
        "core": {"resolution": 0.02, "distance": 2},
        "slab90": {"resolution": 0.05, "distance": 2},
        "Oxide": {"resolution": 0.3, "distance": 2},
    }

    ymin = -30  # make sure all objects cross the line
    ymax = 30  # make sure all objects cross the line
    xcoord = lines_x[0]
    xsection_bounds = [[xcoord, ymin], [xcoord, ymax]]

    mesh_x1 = gf.simulation.gmsh.uz_xsection_mesh(
        c,
        xsection_bounds,
        filtered_layerstack,
        resolutions=resolutions,
        background_tag="Oxide",
        background_padding=(
            2.0,
            2.0,
            2.0,
            2.0,
        ),  # how much of backgorund tag to add to each side of the simulatoin
        filename="mesh_x1.msh",
    )

    ymin = -30  # make sure all objects cross the line
    ymax = 30  # make sure all objects cross the line
    xcoord = lines_x[25]
    xsection_bounds = [[xcoord, ymin], [xcoord, ymax]]

    mesh_x2 = gf.simulation.gmsh.uz_xsection_mesh(
        c,
        xsection_bounds,
        filtered_layerstack,
        resolutions=resolutions,
        background_tag="Oxide",
        background_padding=(
            2.0,
            2.0,
            2.0,
            2.0,
        ),  # how much of backgorund tag to add to each side of the simulatoin
        filename="mesh_x2.msh",
    )
