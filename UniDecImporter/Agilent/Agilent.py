"""High-level Agilent MassHunter ``.d`` importer."""

from __future__ import annotations

import numpy as np

from ..Importer import Importer
from .MZFILE import MZFile


class AgilentImporter(Importer):
    """Read Agilent data through a locally licensed MassHunter Data Access SDK."""

    def __init__(self, path, **kwargs):
        """Open an Agilent ``.d`` directory through MassHunter."""
        super().__init__(path, **kwargs)
        self.msrun = MZFile(path)
        self.init_scans()
        self.cdms_support = False
        self.imms_support = False
        self.chrom_support = True

    def init_scans(self):
        """Load scan identifiers, retention times, and MS levels."""
        info = self.msrun.scan_info()
        self.scans = np.asarray([item[2] + 1 for item in info], dtype=int)
        self.times = np.asarray([item[0] for item in info], dtype=float)
        self.levels = np.asarray([
            2 if "MS2" in item else 1 for item in info
        ], dtype=int)
        self.scan_range = [int(self.scans.min()), int(self.scans.max())]

    def get_all_scans(self, threshold=-1):
        """Load all scans, optionally filtering by intensity."""
        self.data = [self.get_single_scan(scan, threshold=threshold) for scan in self.scans]
        return self.data

    def get_single_scan(self, scan, threshold=-1):
        """Return one scan as a filtered m/z-intensity array."""
        data = np.asarray(self.msrun.scan(int(scan) - 1), dtype=float)
        if data.size == 0:
            return np.empty((0, 2))
        data = data.reshape((-1, 2))
        data = data[data[:, 0] > 10]
        if threshold >= 0:
            data = data[data[:, 1] > threshold]
        return data

    def get_avg_scan(self, scan_range=None, time_range=None, mzbins=None):
        """Merge selected scans, optionally using a fixed m/z bin width."""
        if mzbins is None:
            return self.avg_fast(scan_range, time_range)
        from ..ImportTools import merge_spectra

        scan_range = self.scan_range_from_inputs(scan_range, time_range)
        spectra = [self.get_single_scan(scan) for scan in self.scans
                   if scan_range[0] <= scan <= scan_range[1]]
        return merge_spectra(spectra, mzbins=mzbins)

    def get_tic(self):
        """Return the total-ion chromatogram."""
        return np.asarray(self.msrun.xic(filter="Full"), dtype=float)

    def get_eic(self, mass, mz_tol, rt_range=None):
        """Return an extracted-ion chromatogram using an absolute m/z tolerance."""
        start_time, stop_time = (None, None) if rt_range is None else rt_range
        return np.asarray(self.msrun.xic(
            start_time=start_time,
            stop_time=stop_time,
            start_mz=mass - mz_tol,
            stop_mz=mass + mz_tol,
        ), dtype=float)

    def get_polarity(self, scan=1):
        """Return the polarity recorded for *scan*."""
        item = self.msrun.scan_info()[self.get_scan_index(scan)]
        text = str(item[-1]).lower()
        if "positive" in text or "+" in text:
            return "Positive"
        if "negative" in text or "-" in text:
            return "Negative"
        return None

    def close(self):
        """Close the MassHunter data handle."""
        if getattr(self, "msrun", None) is not None:
            self.msrun.close()
            self.msrun = None

