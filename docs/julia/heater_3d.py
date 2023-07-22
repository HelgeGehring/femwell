import gdsfactory as gf
from gdsfactory.generic_tech import LAYER
from gdsfactory.pdk import LayerStack, get_layer_stack
from gdsfactory.simulation.gmsh.xyz_mesh import xyz_mesh

# Choose some component
c = gf.component.Component()
waveguide = c << gf.get_component(gf.components.straight_heater_metal(length=40))
e1 = c << gf.components.straight(1, cross_section="metal3")
e1.connect(e1["e1"], waveguide["l_e1"])
e2 = c << gf.components.straight(1, cross_section="metal3")
e2.connect(e2["e1"], waveguide["r_e3"])
c.add_port("e1", port=e1["e2"])
c.add_port("e2", port=e2["e2"])
c.show()

# Add wafer / vacuum (could be automated)
wafer = c << gf.components.bbox(bbox=waveguide.bbox, layer=LAYER.WAFER)

# Generate a new component and layerstack with new logical layers
layerstack = get_layer_stack()
c = layerstack.get_component_with_net_layers(
    c,
    portnames=["e1", "e2"],
    delimiter="#",
)

# FIXME: .filtered returns all layers
# filtered_layerstack = layerstack.filtered_from_layerspec(layerspecs=c.get_layers())
filtered_layerstack = LayerStack(
    layers={
        k: layerstack.layers[k]
        for k in (
            # "via1",
            "box",
            "clad",
            # "metal2",
            "metal3#e1",
            "heater",
            "via2",
            "core",
            "metal3#e2",
            "metal3",
            "metal2",
            "via1"
            # "metal3",
            # "via_contact",
            # "metal1"
        )
    }
)
print(layerstack.layers.keys())

filtered_layerstack.layers["clad"].thickness += sum(
    filtered_layerstack.layers[name].thickness for name in ["via2", "metal3"]
)

resolutions = {
    "core": {"resolution": 0.3},
    "via2": {"resolution": 0.3},
    "via1": {"resolution": 0.3},
    "heater": {"resolution": 0.3},
}
geometry = xyz_mesh(
    component=c,
    layerstack=filtered_layerstack,
    resolutions=resolutions,
    filename="mesh.msh",
    default_characteristic_length=2,
    verbosity=5,
    global_scaling=1e-6,
)
