"""Package-specific exceptions."""


class UniDecImporterError(Exception):
    """Base error for this package."""


class UnsupportedFormatError(UniDecImporterError, OSError):
    """Raised when no reader is registered for a path."""


class VendorReaderUnavailableError(UniDecImporterError):
    """Raised when a vendor reader cannot run on the current platform."""

