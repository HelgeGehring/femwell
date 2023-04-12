"""
Adapted from F. Laporte at https://github.com/flaport/meow/blob/main/meow/eme/common.py
and references.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

from femwell.mode_solver import (
    calculate_hfield,
    calculate_overlap,
    calculate_scalar_product,
    plot_mode,
)

sys.path.insert(0, "./mesh")
import sax

from femwell.mesh.slice import slice_component_xbounds

try:
    import klujax
except ImportError:
    klujax = None


def compute_interface_s_matrix(
    mode_a,  # lams, basis, xs
    mode_b,
):
    lams_a, basis_a, xs_a = mode_a
    lams_b, basis_b, xs_b = mode_b

    products_ab = np.zeros((len(lams_a), len(lams_b)), dtype=complex)
    products_ba = np.zeros((len(lams_b), len(lams_a)), dtype=complex)
    for i in range(len(lams_a)):
        for j in range(len(lams_b)):
            E_i = xs_a[i]
            E_j = xs_b[j]
            H_i = calculate_hfield(basis_a, E_i, lams_a[i] * (2 * np.pi / 1.55))
            H_j = calculate_hfield(basis_b, E_j, lams_b[j] * (2 * np.pi / 1.55))
            products_ab[i, j] = np.abs(calculate_scalar_product(basis_a, E_i, basis_b, H_j))
            products_ba[j, i] = np.abs(calculate_scalar_product(basis_b, E_j, basis_a, H_i))

    T_ab = 2 * np.linalg.inv(products_ab + products_ba.T)
    R_ab = 0.5 * (products_ba.T - products_ab) @ T_ab

    T_ba = 2 * np.linalg.inv(products_ba + products_ab.T)
    R_ba = 0.5 * (products_ab.T - products_ba) @ T_ba

    S = np.concatenate(
        [
            np.concatenate([R_ab, T_ba], 1),
            np.concatenate([T_ab, R_ba], 1),
        ],
        0,
    )

    # create port map
    in_ports = [f"left@{i}" for i in range(len(mode_a))]
    out_ports = [f"right@{i}" for i in range(len(mode_b))]
    port_map = {p: i for i, p in enumerate(in_ports + out_ports)}
    return S, port_map
    # return T_ab, R_ab, T_ba, R_ba


def compute_propagation_s_matrix(modes, length, wavelength):
    lams, basis, xs = modes  # lams is neff
    betas = lams * 2 * np.pi / wavelength
    # S = np.diag(np.exp(-1 * 2j * np.pi * np.abs(betas) * length))
    s_dict = {(f"left@{i}", f"right@{i}"): np.exp(beta * length) for i, beta in enumerate(betas)}
    s_dict = {**s_dict, **{(p2, p1): v for (p1, p2), v in s_dict.items()}}
    return s_dict


def _get_netlist(propagations, interfaces):
    """get the netlist of a stack of `Modes`"""

    instances = {
        **{k: _load_constant_model(S) for k, S in propagations.items()},
        **{k: _load_constant_model(S) for k, S in interfaces.items()},
    }

    connections = {}
    for i in range(len(interfaces)):
        connections[f"p_{i},right"] = f"i_{i}_{i+1},left"
        connections[f"i_{i}_{i+1},right"] = f"p_{i+1},left"

    ports = {
        f"left": f"p_0,left",
        f"right": f"p_{len(propagations)-1},right",
    }

    return {"instances": instances, "connections": connections, "ports": ports}


def _validate_sax_backend(sax_backend):
    if sax_backend is None:
        sax_backend = "klu" if klujax is not None else "default"

    if sax_backend not in ["default", "klu"]:
        raise ValueError(
            f"Invalid SAX Backend. Got: {sax_backend!r}. Should be 'default' or 'klu'."
        )
    return sax_backend


def _load_constant_model(value):
    def model():
        return value

    return model


def compute_total_S_matrix(
    meshnames,
    mesh_info_dict,
    lengths,
    wavelength,
    num_modes,
    sax_backend=None,
):
    """
    Given a list of N meshes corresponding to sections lengths, computes the overall S-matrix of the resulting modes.
    Uses a propagation matrix along each length, and an interface matrix between each segment.

        length[0]              length[1]                                   length[N]
    <---------------><----------------------------->            <----------------------------->
    |---meshes[0]---|----------meshes[1]-----------|---.....----|----------meshes[N]-----------|
    ^               ^                              ^            ^                              ^
    In              0-1                            1-2        (N-1)-N                          Out
                  interface                     interface    interface

        Args:
            meshes: list of gmsh mesh objects, with physicals named like keys of indices
            indices_dict: Dict[physical names: refractive index]
            lengths: list of float propagation lengths between interfaces
            wavelength: float
            num_modes: int, number of modes to compute for each eigensolve
    """
    modes = []
    propagations = {}
    for i, meshname in enumerate(meshnames):
        # with "." as tmpdirname: # tempfile.TemporaryDirectory()
        mesh = Mesh.load(meshname)
        basis0 = Basis(mesh, ElementTriP0(), intorder=4)
        epsilon = basis0.zeros()
        for name, refractive_index in indices_dict.items():
            try:
                epsilon[basis0.get_dofs(elements=name)] = refractive_index**2
            except ValueError:
                pass
        lams, basis, xs = compute_modes(
            basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=num_modes
        )
        plot_mode(basis, np.real(xs[0]))
        plt.show()
        plot_mode(basis, np.imag(xs[0]))
        plt.show()
        modes.append((lams, basis, xs))
        # Create SAX model for propagation
        propagations[f"p_{i}"] = compute_propagation_s_matrix(
            (lams, basis, xs), lengths[i], wavelength
        )
    interfaces = {
        f"i_{i}_{i + 1}": compute_interface_s_matrix(modes[i], modes[i + 1])
        for i in range(len(meshnames) - 1)
    }

    net = _get_netlist(propagations, interfaces)
    mode_names = [f"{i}" for i in range(num_modes)]

    sax_backend = _validate_sax_backend(sax_backend)
    _circuit, _ = sax.circuit(netlist=net, backend=sax_backend, modes=mode_names)

    # maybe in the future we should return the sax model and not the S-matrix?
    sdict = sax.sdict(_circuit())
    sdict = {k: sdict[k] for k in sorted(sdict)}
    return sax.sdense(sdict)


if __name__ == "__main__":
    import tempfile

    # from femwell.waveguide import mesh_waveguide
    import gdsfactory as gf
    import gdsfactory.simulation.gmsh as gfmesh
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm
    from skfem import Basis, ElementTriP0, Mesh
    from tqdm.auto import tqdm

    from femwell.mode_solver import (
        calculate_hfield,
        calculate_overlap,
        compute_modes,
        plot_mode,
    )

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

    # taper = c.add_ref(
    #     gf.components.taper(length = 10.0,
    #                             width1 = 0.5,
    #                             width2 = 2,
    #                             cross_section = "rib",
    #             )
    # )
    # straight = c.add_ref(gf.components.straight(10, width = 2, cross_section="rib"))
    # straight.connect("o1", taper.ports["o2"])
    straight2 = c.add_ref(gf.components.straight(10, width=2, cross_section="strip"))
    # straight2.connect("o1", straight.ports["o2"])
    straight = c.add_ref(gf.components.straight(10, width=2, cross_section="rib"))
    straight.connect("o1", straight2.ports["o2"])
    # taper = c.add_ref(
    #     gf.components.taper(length = 10.0,
    #                             width1 = 0.5,
    #                             width2 = 2,
    #                             cross_section = "rib",
    #             )
    # )
    # taper.connect("o2", straight.ports["o2"])
    # bend = c.add_ref(
    #     gf.components.bend_euler(angle = 20.0,
    #                             p=0.5,
    #                             width=0.5,
    #                             cross_section = "strip",
    #             ).move([12, -5])
    # )

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
    nm = 1e-3
    # Compute the x-coordinates where the cross-section changes
    lines_x = slice_component_xbounds(c, filtered_layerstack, mesh_step=1000 * nm)
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

    ymin = -6  # make sure all objects cross the line
    ymax = 6  # make sure all objects cross the line

    meshes = []
    lengths = []
    meshnames = []
    for x1, x2 in zip(lines_x[:-1], lines_x[1:]):
        x = (x1 + x2) / 2
        lengths.append(x2 - x1)
        xsection_bounds = [[x, ymin], [x, ymax]]
        meshes.append(
            gf.simulation.gmsh.uz_xsection_mesh(
                c,
                xsection_bounds,
                filtered_layerstack,
                resolutions=resolutions,
                background_tag="clad",
                background_padding=(
                    2.0,
                    2.0,
                    2.0,
                    2.0,
                ),  # how much of backgorund tag to add to each side of the simulatoin
                filename=f"mesh_x_{x1}.msh",
            )
        )
    meshnames.append(f"mesh_x_{x1}.msh")
    num_modes = 16

    indices_dict = {"core": 3.45, "slab90": 3.45, "clad": 1.44}

    S = compute_total_S_matrix(
        meshnames, indices_dict, lengths, wavelength=1.55, num_modes=num_modes
    )
    print(S)
