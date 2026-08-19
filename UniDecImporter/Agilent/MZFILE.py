"""Small pythonnet wrapper around the Agilent MassHunter Data Access SDK.

The interface originated in UniDec's multiplierz-derived Agilent reader and is kept
private to the high-level :class:`AgilentImporter`.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

try:
    import clr
except ImportError as error:  # pragma: no cover - exercised on Windows without the extra
    raise ImportError("Agilent support requires the 'vendor' or 'agilent' extra (pythonnet)") from error


SDK_DIR = Path(os.environ.get("AGILENT_DA_SDK_DIR", Path(__file__).parent))
ASSEMBLIES = ("MassSpecDataReader", "BaseCommon", "BaseDataAccess", "MIDAC", "agtsampleinforw")

for assembly in ASSEMBLIES:
    assembly_path = SDK_DIR / f"{assembly}.dll"
    try:
        clr.AddReference(str(assembly_path) if assembly_path.is_file() else assembly)
    except Exception as error:  # pragma: no cover - depends on the licensed SDK
        raise OSError(
            "Unable to load the Agilent MassHunter Data Access SDK. Set "
            "AGILENT_DA_SDK_DIR to its assembly directory."
        ) from error

clr.AddReference("System.Collections")

from Agilent.MassSpectrometry.DataAnalysis import (  # type: ignore[import-not-found]  # noqa: E402
    BDAChromFilter,
    ChromType,
    DesiredMSStorageType,
    IBDAChromFilter,
    IMsdrChargeStateAssignmentFilter,
    IMsdrDataReader,
    MassSpecDataReader,
    MinMaxRange,
    MSLevel,
    MsdrChargeStateAssignmentFilter,
    MsdrPeakFilter,
)


for interface, implementation in (
    (IMsdrDataReader, MassSpecDataReader),
    (IBDAChromFilter, BDAChromFilter),
    (IMsdrChargeStateAssignmentFilter, MsdrChargeStateAssignmentFilter),
):
    for name in dir(interface):
        try:
            setattr(implementation, name, getattr(interface, name))
        except Exception:
            pass


class MZFile:
    """Read one Agilent ``.d`` directory through the vendor SDK."""

    def __init__(self, datafile, **kwargs):
        del kwargs
        self.file_type = ".d"
        self.data_file = os.fspath(datafile)
        self._scan_info = None
        self.tic_object = None
        self.no_filter = MsdrPeakFilter()
        self.source = MassSpecDataReader()
        try:
            opened = self.source.OpenDataFile(self.source, self.data_file)
        except Exception as error:
            raise OSError(f"error opening Agilent data: {self.data_file}") from error
        if not opened:
            raise OSError(f"error opening Agilent data: {self.data_file}")

    def close(self):
        if self.source is not None:
            self.source.CloseDataFile(self.source)
            self.source = None

    def scan_range(self):
        return 0, int(self.source.FileInformation.MSScanFileInformation.TotalScansPresent)

    def scan_info(self, start_time=None, stop_time=None, start_mz=None, stop_mz=None):
        if self._scan_info is None:
            self._scan_info = []
            for index in range(self.scan_range()[1]):
                record = self.source.GetScanRecord(self.source, index)
                self._scan_info.append((
                    float(record.RetentionTime),
                    float(record.MZOfInterest),
                    index,
                    f"MS{int(record.MSLevel)}",
                    str(record.IonPolarity),
                ))
        result = self._scan_info
        if start_time is not None:
            result = [item for item in result if item[0] > start_time]
        if stop_time is not None:
            result = [item for item in result if item[0] < stop_time]
        if start_mz is not None:
            result = [item for item in result if item[1] > start_mz]
        if stop_mz is not None:
            result = [item for item in result if item[1] < stop_mz]
        return result

    def headers(self):
        return self.scan_info()

    def scan(self, index, mode=None):
        modes = {
            None: DesiredMSStorageType.PeakElseProfile,
            "peakelseprofile": DesiredMSStorageType.PeakElseProfile,
            "profileelsepeak": DesiredMSStorageType.ProfileElsePeak,
            "profile": DesiredMSStorageType.Profile,
            "peak": DesiredMSStorageType.Peak,
            "centroid": DesiredMSStorageType.Peak,
        }
        requested = None if mode is None else mode.lower()
        if requested not in modes:
            raise ValueError(f"unsupported Agilent spectrum mode: {mode}")
        spectrum = self.source.GetSpectrum(
            self.source, int(index), self.no_filter, self.no_filter, modes[requested]
        )
        return list(zip(spectrum.XArray, spectrum.YArray))

    def xic(self, start_time=None, stop_time=None, start_mz=None, stop_mz=None,
            filter=None, UV=False):
        if filter and filter.strip().lower() not in {"full", "full ms"}:
            raise ValueError("only full-MS chromatograms are supported")
        start_time = 0 if start_time is None else start_time
        stop_time = 999999 if stop_time is None else stop_time
        start_mz = 0 if start_mz is None else start_mz
        stop_mz = 999999 if stop_mz is None else stop_mz

        chrom_filter = BDAChromFilter()
        chrom_filter.set_MSLevelFilter(chrom_filter, MSLevel.MS)
        chrom_filter.set_ChromatogramType(
            chrom_filter, ChromType.ExtractedWavelength if UV else ChromType.ExtractedIon
        )
        chrom_filter.set_SingleChromatogramForAllMasses(chrom_filter, True)

        mz_range = MinMaxRange()
        mz_range.set_Min(float(start_mz))
        mz_range.set_Max(float(stop_mz))
        chrom_filter.set_IncludeMassRanges(chrom_filter, (mz_range,))

        rt_range = MinMaxRange()
        rt_range.set_Min(float(start_time))
        rt_range.set_Max(float(stop_time))
        chrom_filter.set_ScanRange(chrom_filter, rt_range)

        chromatogram = self.source.GetChromatogram(self.source, chrom_filter).Get(0)
        return list(zip(chromatogram.XArray, chromatogram.YArray))

    def uv_trace(self):
        devices = self.source.GetNonmsDevices()
        if not devices.Length:
            raise OSError("no non-MS devices were available")
        return self.source.GetTWC(devices[0])
