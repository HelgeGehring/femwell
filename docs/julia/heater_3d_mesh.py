import gdsfactory as gf
from gdsfactory.generic_tech import LAYER
from gdsfactory.pdk import LayerStack, get_layer_stack
from gplugins.common.utils.get_component_with_net_layers import (
    get_component_with_net_layers,
)
from gplugins.gmsh.get_mesh import get_mesh

# Choose some component
c = gf.component.Component()
waveguide = c << gf.get_component(gf.components.straight_heater_metal(length=40))
e1 = c << gf.components.straight(1, cross_section="xs_m3")
e1.connect(e1["e1"], waveguide["l_e1"])
e2 = c << gf.components.straight(1, cross_section="xs_m3")
e2.connect(e2["e1"], waveguide["r_e3"])
c.add_port("e1", port=e1["e2"])
c.add_port("e2", port=e2["e2"])
c.show()

# Add wafer / vacuum (could be automated)
wafer = c << gf.components.bbox(bbox=waveguide.bbox, layer=LAYER.WAFER)

# Generate a new component and layerstack with new logical layers
layerstack = get_layer_stack().filtered(
    ("box", "clad", "heater", "via2", "core", "metal3", "metal2", "via1")
)
c = get_component_with_net_layers(
    c,
    layerstack,
    port_names=["e1", "e2"],
    delimiter="#",
)

layerstack.layers["clad"].thickness += sum(
    layerstack.layers[name].thickness for name in ["via2", "metal3"]
)

resolutions = {
    "core": {"resolution": 0.3},
    "via2": {"resolution": 0.2},
    "via1": {"resolution": 0.2},
    "heater": {"resolution": 0.5},
}
geometry = get_mesh(
    type="3D",
    component=c,
    layer_stack=layerstack,
    resolutions=resolutions,
    filename="mesh.msh",
    default_characteristic_length=2,
    # global_3D_algorithm=10,
    verbosity=5,
    global_scaling=1e-6,
)
