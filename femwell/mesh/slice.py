from gdsfactory.simulation.gmsh import fuse_component_layer, order_layerstack
from collections import OrderedDict


def get_layer_bounds(polygon):
    """Process a polygon to identify where it does not change along propagation direction.
    
        Propagation direction is "x" in software.
    """
    # for vertex in polygon.bounds.exterior:
    return True


def slice_component():
    """Returns list of z-coordinates where cross-section is to be taken."""
    # for layer in 
    return True


def overlap_mesh():
    """Returns a mesh conditioned on the shapes from two different cross-sections."""
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

    # Reorder polygons according to meshorder
    layer_order = order_layerstack(filtered_layerstack)
    ordered_layers = [value for value in layer_order if value in set(layer_order)]
    shapes = OrderedDict()
    for layer in ordered_layers:
        shapes[layer] = layer_polygons_dict[layer]

    # Compute the x-coordinates where the cross-section changes

    polygon = gf.simulation.gmsh.get_xsection_bound_polygons(component=c, 
                                                    xsection_bounds=[[5, -200], [5, 200]], 
                                                    layer_stack=filtered_layerstack
                                                )

    print(polygon)

