"""Package-specific exceptions."""


class MassSpecImporterError(Exception):
    """Base error for this package."""


class UnsupportedFormatError(MassSpecImporterError, OSError):
    """Raised when no reader is registered for a path."""


class VendorReaderUnavailableError(MassSpecImporterError):
    """Raised when a vendor reader cannot run on the current platform."""

