"""Standalone readers for mass-spectrometry data."""

from ._version import __version__
from .exceptions import UnsupportedFormatError, VendorReaderUnavailableError
from .ImporterFactory import ImporterFactory, get_polarity, recognized_types

def get_importer(file_path, **kwargs):
    """Return the format-specific reader selected for *file_path*."""
    return ImporterFactory.create_importer(file_path, **kwargs)


__all__ = [
    "ImporterFactory",
    "UnsupportedFormatError",
    "VendorReaderUnavailableError",
    "__version__",
    "get_importer",
    "get_polarity",
    "recognized_types",
]

