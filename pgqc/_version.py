from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pgqc")
except PackageNotFoundError:
    __version__ = "uninstalled"
