from femwell.mesh.meshtracker import MeshTracker
from femwell.mesh.mesh import mesh_from_OrderedDict, break_line
from femwell.mesh.slice import process_component, get_vertices, get_polygon_x_bounds, get_unique_x_bounds, get_mode_regions, slice_component, overlap_mesh
__all__ = [
    "MeshTracker",
    "mesh_from_OrderedDict",
    "break_line",
    "process_component", 
    "get_vertices", 
    "get_polygon_x_bounds", 
    "get_unique_x_bounds", 
    "get_mode_regions", 
    "slice_component", 
    "overlap_mesh",
]
__version__ = "0.0.1"
