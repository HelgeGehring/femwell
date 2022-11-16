.. include:: ../../README.rst

.. toctree::
    thermal.rst
    maxwell.rst



.. plot::
    import tempfile

    from collections import OrderedDict

    import matplotlib.pyplot as plt
    import numpy as np
    import shapely.geometry

    from skfem import Mesh, Basis, ElementTriP0, ElementTriP1, Functional

    from femwell.mode_solver import compute_modes, plot_mode
    from femwell.mesh import mesh_from_OrderedDict

    with tempfile.TemporaryDirectory() as tmpdirname:
        core = shapely.geometry.box(-.1, -.15, .1, .15)

        polygons = OrderedDict(
            core=core,
            clad=core.buffer(15, resolution=4)
        )

        resolutions = dict(
            core={"resolution": .01, "distance": .1}
        )

        mesh_from_OrderedDict(polygons, resolutions, filename='mesh.msh', default_resolution_max=10)
        mesh_from_OrderedDict(polygons, resolutions, filename=tmpdirname + '/mesh.msh', default_resolution_max=10)
        mesh = Mesh.load(tmpdirname + '/mesh.msh')

    basis0 = Basis(mesh, ElementTriP0(), intorder=4)
    epsilon = basis0.zeros().astype(complex)
    epsilon[basis0.get_dofs(elements='core')] = 1.9963 ** 2
    epsilon[basis0.get_dofs(elements='clad')] = 1.444 ** 2
    # basis0.plot(np.real(epsilon), colorbar=True).show()

    lams, basis, xs = compute_modes(basis0, epsilon, wavelength=1.55, mu_r=1, num_modes=1)
    basis.plot(np.real(xs[0]), colorbar=True).show()