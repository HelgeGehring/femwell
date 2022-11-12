import gdsfactory as gf
from gdsfactory.simulation.gmsh import fuse_component_layer, order_layerstack
from gdsfactory.geometry.boolean import boolean
from collections import OrderedDict
import numpy as np

from shapely.geometry import MultiPoint, Polygon, LineString
from shapely.affinity import translate

nm = 1E-3

def process_component(component, layerstack):
    """Process component polygons to:
        * eliminate precision issues on vertices
        * fuse polygons into the smallest set of polygons

    Returns
        layer_polygons_dict: dict containing layername as key, and simplest MultiPolygon object as entry
    """
    layer_dict = layerstack.to_dict()
    layer_polygons_dict = {}
    for layername in layer_dict.keys():
        layer_polygons_dict[layername] = fuse_component_layer(
            c, layername, layer_dict[layername]
        )
    return layer_polygons_dict

def get_vertices(polygon):
    """Return all polygon vertices (interior and exterior)"""
    vertices = []
    for polygon_hole in list(polygon.interiors):
        for vertex in MultiPoint(polygon_hole.coords).geoms:
            vertices.append(vertex)
    # Parse boundary
    for vertex in MultiPoint(
        polygon.exterior.coords
    ).geoms:
        vertices.append(vertex)

    return vertices

def get_polygon_x_bounds(polygon):
    """Return x_bounds of polygon vertices.
    
        Propagation direction is "x" in component ("z" in mode solver)
    """
    xs = [p.x for p in get_vertices(polygon)]
    return xs

def get_unique_x_bounds(layer_polygons_dict):
    """Return unique x_bounds across all polygon vertices of a layer_polygons_dict
    
        Propagation direction is "x" in component ("z" in mode solver)
    """
    xs = []
    for polygons in layer_polygons_dict.values():
        for polygon in polygons.geoms if hasattr(polygons, "geoms") else [polygons]:
            xs.extend(get_polygon_x_bounds((polygon)))
    return np.sort(np.unique(np.array(xs)))

def get_mode_regions(component, layerstack, tol=1E-6):
    """Return interesting x_bounds of polygon vertices.
    
        Propagation direction is "x" in component ("z" in mode solver)

        Args:
            component: gdsfactory Component to process
            layerstack: gdsfactory LayerStack to process
        Returns:
            x_changing_bounds: x-bounds where the component needs to be meshed
            x_not_changing_bounds: x-bounds where there the structure is unchanging (free propagation)
    """
    bbox = gf.components.bbox(bbox=component.bbox).get_polygons()[0][:,1]
    ymin = np.min(bbox)
    ymax = np.max(bbox)
    layer_polygons_dict =  process_component(component, layerstack)
    x_bounds = get_unique_x_bounds(layer_polygons_dict)

    x_changing_bounds = []
    x_not_changing_bounds = []
    for x1, x2 in [
            (x_bounds[i] + tol, x_bounds[i + 1] - tol) for i in range(0, len(x_bounds) - 1)
        ]:
        found_different = False
        for layername, polygons in layer_polygons_dict.items():
            if found_different:
                continue
            else:
                line_x1 = LineString([[x1, ymin], [x1, ymax]])
                line_x2 = LineString([[x2, ymin], [x2, ymax]])

                xsection_x1 = polygons.intersection(line_x1)
                xsection_x2 = polygons.intersection(line_x2)

                if not xsection_x1.equals(translate(xsection_x2, xoff=x1-x2)):
                    found_different = True
                    x_changing_bounds.append([x1-tol, x2+tol])

        if not found_different:
            x_not_changing_bounds.append([x1-tol, x2+tol])

    return x_changing_bounds, x_not_changing_bounds


def slice_component(component, layerstack, mesh_step=100*nm):
    """Returns minimal list of x-coordinates where cross-section is to be taken."""
    # Process polygon to extract regions of change and free propagation
    x_changing_bounds, x_not_changing_bounds = get_mode_regions(component, layerstack, tol=1E-6)

    # Where geometry is changing, mesh according to mesh_step
    x_coordinates = []
    for x1, x2 in x_changing_bounds:
        x_coordinates.append(np.arange(x1, x2, mesh_step))
    # Where not changing, just return the bounds
    for x1, x2 in x_not_changing_bounds:
        x_coordinates.append(np.array([x1]))

    # Return sorted bounds
    return np.sort(np.concatenate(x_coordinates).ravel())


def overlap_mesh(component, layerstack, mesh_step=100*nm):
    return True
#     """Returns a mesh conditioned on the shapes from N different cross-sections."""
#     # Get cross-sectional profiles from each x-coordinate mesh point
#     x_coords = slice_component(component, layerstack, mesh_step=100*nm)
#     shapes = {}
#     for x_coord in x_coords:
#         shapes[x] = 
#     return True


if __name__ == "__main__":

    import gdsfactory as gf
    import gdsfactory.simulation.gmsh as gfmesh

    c = gf.Component()

    taper = c.add_ref(
        gf.components.taper(length = 10.0,
                                width1 = 0.5,
                                width2 = 2,
                                cross_section = "rib",
                )
    )
    straight = c.add_ref(gf.components.straight(10, width = 2, cross_section="rib"))
    straight.connect("o1", taper.ports["o2"])
    straight2 = c.add_ref(gf.components.straight(10, width = 2, cross_section="strip"))
    straight2.connect("o1", straight.ports["o2"])
    straight = c.add_ref(gf.components.straight(10, width = 2, cross_section="rib"))
    straight.connect("o1", straight2.ports["o2"])
    taper = c.add_ref(
        gf.components.taper(length = 10.0,
                                width1 = 0.5,
                                width2 = 2,
                                cross_section = "rib",
                )
    )
    taper.connect("o2", straight.ports["o2"])
    bend = c.add_ref(
        gf.components.bend_euler(angle = 20.0,
                                p=0.5,
                                width=0.5,
                                cross_section = "strip",
                ).move([12, -5])
    )

    from gdsfactory.tech import get_layer_stack_generic, LayerStack
    import numpy as np

    filtered_layerstack = LayerStack(
        layers={
            k: get_layer_stack_generic().layers[k]
            for k in (
                "slab90",
                "core",
            )
        }
    )

    # Get unique cross-sections

    # Reorder polygons according to meshorder
    # layer_order = order_layerstack(filtered_layerstack)
    # ordered_layers = [value for value in layer_order if value in set(layer_order)]
    # shapes = OrderedDict()
    # for layer in ordered_layers:
    #     shapes[layer] = layer_polygons_dict[layer]

    # Compute the x-coordinates where the cross-section changes
    lines = slice_component(c, filtered_layerstack)
    for x in lines:
        P = gf.Path([[x, -20], [x, 20]])
        X = gf.CrossSection(width=0.001, layer=(99,0))
        line = gf.path.extrude(P, X)
        c << line

    c.show()


    # polygons_dict = gf.simulation.gmsh.get_xsection_bound_polygons(component=c, 
    #                                                 xsection_bounds=[[5, -200], [5, 200]], 
    #                                                 layer_stack=filtered_layerstack
    #                                             )

    # for name, polygons in polygons_dict.items():
    #     for polygon in polygons:
    #         print(name, get_layer_bounds(polygon))

