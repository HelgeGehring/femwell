try:
    from importlib import metadata
except ImportError:  # for Python<3.8
    import importlib_metadata as metadata

try:
    __version__ = metadata.version(__package__ or __name__)
except metadata.PackageNotFoundError:
    __version__ = "git"
