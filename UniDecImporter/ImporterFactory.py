"""Select an importer from a mass-spectrometry data path."""

from __future__ import annotations

import os
import platform
from pathlib import Path

from . import SingleScanImporter as SSI
from .I2MS.I2MS import I2MSImporter
from .MZML.mzML import MZMLImporter
from .MZXML.mzXML import MZXMLImporter
from .exceptions import UnsupportedFormatError, VendorReaderUnavailableError


OPEN_FORMATS = (".mzxml", ".mzml", ".mzml.gz", ".gz", ".txt", ".dat", ".csv",
                ".npz", ".i2ms", ".dmt", ".bin")
VENDOR_FORMATS = (".raw", ".d")
recognized_types = list(VENDOR_FORMATS + OPEN_FORMATS)


def _is_windows_x64() -> bool:
    """Return whether the current platform can load the supported vendor SDKs."""
    return platform.system() == "Windows" and platform.machine().lower() in {
        "amd64", "x86_64", "x64"
    }


def _vendor_platform_error(vendor: str) -> VendorReaderUnavailableError:
    """Build a platform error naming the unavailable *vendor* reader."""
    return VendorReaderUnavailableError(
        f"{vendor} data require 64-bit x86 Windows and the bundled vendor runtime; "
        f"current platform is {platform.system()} {platform.machine()}"
    )


class ImporterFactory:
    """Create the reader appropriate for a path's extension and shape."""

    def __init__(self):
        """Initialize the factory's list of recognized file extensions."""
        self.recognized_file_types = list(recognized_types)

    @staticmethod
    def create_importer(file_path, **kwargs):
        """Create and return the reader appropriate for *file_path*.

        Parameters
        ----------
        file_path : str or os.PathLike
            Mass-spectrometry file, or a Waters ``.raw``/Agilent ``.d`` directory.
        **kwargs
            Options forwarded to the selected reader.

        Raises
        ------
        UnsupportedFormatError
            If the path extension is not recognized.
        VendorReaderUnavailableError
            If a proprietary reader cannot run on the current platform.
        """
        path = Path(file_path)
        name = path.name.lower()
        ending = ".mzml.gz" if name.endswith(".mzml.gz") else path.suffix.lower()

        if ending == ".raw":
            vendor = "Waters" if path.is_dir() else "Thermo"
            if not _is_windows_x64():
                raise _vendor_platform_error(vendor)
            try:
                if path.is_dir():
                    from .Waters.Waters import WatersDataImporter

                    return WatersDataImporter(str(path), **kwargs)
                from .Thermo.Thermo import ThermoImporter

                return ThermoImporter(str(path), **kwargs)
            except (ImportError, OSError) as error:
                raise VendorReaderUnavailableError(
                    f"The {vendor} reader could not load its vendor runtime: {error}"
                ) from error

        if ending == ".d":
            if not _is_windows_x64():
                raise _vendor_platform_error("Agilent")
            try:
                from .Agilent.Agilent import AgilentImporter

                return AgilentImporter(str(path), **kwargs)
            except (ImportError, OSError) as error:
                raise VendorReaderUnavailableError(
                    f"The Agilent reader could not load its vendor runtime: {error}"
                ) from error
        if ending == ".mzxml":
            return MZXMLImporter(str(path), **kwargs)
        if ending in {".mzml", ".mzml.gz", ".gz"}:
            return MZMLImporter(str(path), **kwargs)
        if ending in {".txt", ".dat", ".csv", ".npz", ".bin"}:
            return SSI.SingleScanImporter(str(path), **kwargs)
        if ending in {".dmt", ".i2ms"}:
            return I2MSImporter(str(path))
        raise UnsupportedFormatError(f"unsupported file type {ending!r}: {file_path}")


def get_polarity(path):
    """Return the detected polarity, defaulting to ``'Positive'`` on read failure."""
    try:
        with ImporterFactory.create_importer(path) as importer:
            return importer.get_polarity() or "Positive"
    except Exception:
        return "Positive"
