import gdsfactory as gf
from gdsfactory.simulation.gmsh import fuse_component_layer, order_layerstack
from gdsfactory.geometry.boolean import boolean
from collections import OrderedDict
import numpy as np

from shapely.geometry import MultiPoint, Polygon, LineString


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

def get_component_x_bounds(layer_polygons_dict):
    """Return unique x_bounds across all polygon vertices of a layer_polygons_dict
    
        Propagation direction is "x" in component ("z" in mode solver)
    """
    xs = []
    for polygons in layer_polygons_dict.values():
        for polygon in polygons.geoms if hasattr(polygons, "geoms") else [polygons]:
            xs.extend(get_polygon_x_bounds((polygon)))
    return np.sort(np.unique(np.array(xs)))

def get_mode_regions(component, layerstack, ymin = -200, ymax = 200, line_width = 1E-4, line_layer = (99, 0)):
    """Return interesting x_bounds of polygon vertices.
    
        Propagation direction is "x" in component ("z" in mode solver)
    """
    layer_polygons_dict =  process_component(component, layerstack)
    x_bounds = get_component_x_bounds(layer_polygons_dict)

    x_parsed_bounds = []
    for x1, x2 in [
            (x_bounds[i], x_bounds[i + 1]) for i in range(0, len(x_bounds) - 1)
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
                
                if not xsection_x1.equals(xsection_x2):
                    found_different = True
                    x_parsed_bounds.append(x1)

    print(x_parsed_bounds)

def slice_component():
    """Returns list of x-coordinates where cross-section is to be taken."""
    # for layer in 
    return True


def overlap_mesh():
    """Returns a mesh conditioned on the shapes from N different cross-sections."""
    return True


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
    straight = c.add_ref(gf.components.straight(20, width = 2, cross_section="rib"))
    straight.connect("o1", taper.ports["o2"])
    taper = c.add_ref(
        gf.components.taper(length = 10.0,
                                width1 = 0.5,
                                width2 = 2,
                                cross_section = "rib",
                )
    )
    taper.connect("o2", straight.ports["o2"])

    c.show()

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

    # Fuse and cleanup polygons of same layer in case user overlapped them
    layer_dict = filtered_layerstack.to_dict()
    layer_polygons_dict = {}
    for layername in layer_dict.keys():  # filtered_layerdict.items():
        layer_polygons_dict[layername] = fuse_component_layer(
            c, layername, layer_dict[layername]
        )

    # Get unique cross-sections

    # Reorder polygons according to meshorder
    layer_order = order_layerstack(filtered_layerstack)
    ordered_layers = [value for value in layer_order if value in set(layer_order)]
    shapes = OrderedDict()
    for layer in ordered_layers:
        shapes[layer] = layer_polygons_dict[layer]

    # Compute the x-coordinates where the cross-section changes
    get_mode_regions(c, filtered_layerstack, ymin = -200, ymax = 200, line_width = 1E-3, line_layer = (99, 0))

    # polygons_dict = gf.simulation.gmsh.get_xsection_bound_polygons(component=c, 
    #                                                 xsection_bounds=[[5, -200], [5, 200]], 
    #                                                 layer_stack=filtered_layerstack
    #                                             )

    # for name, polygons in polygons_dict.items():
    #     for polygon in polygons:
    #         print(name, get_layer_bounds(polygon))

