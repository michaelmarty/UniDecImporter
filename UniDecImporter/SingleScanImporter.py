"""Reader for spectra stored directly in text, NumPy, or binary arrays."""

import os
import numpy as np
from .Importer import Importer
from .ImportTools import header_test


class SingleScanImporter(Importer):
    """Read a single spectrum from TXT, DAT, CSV, NPZ, or BIN data."""

    def __init__(self, filename, **kwargs):
        """Open *filename* and load its array data."""
        super().__init__(filename, **kwargs)

        self.scans = [1]
        self.times = [0]
        self.scan_range = [1, 1]
        self.centroided = False
        self.polarity = "Positive"
        self.scan_number = 1
        self.injection_time = 1
        self.load_data()
        self.cdms_support = True
        self.imms_support = True
        self.chrom_support = False

    def __len__(self):
        """Return the number of rows in the loaded array."""
        return len(self.data)

    def load_data(self):
        """Load, validate, cache, and return the source array."""
        if self.data is not None:
            return self.data
        ext = self.ext.lower()
        if ext == ".txt" or ext == ".dat":
            self.data = np.loadtxt(self._file_path, skiprows=header_test(self._file_path))
        elif ext == ".csv":
            self.data = np.loadtxt(self._file_path, delimiter=",", skiprows=header_test(self._file_path),
                                   usecols=(0, 1))
        elif ext == '.npz':
            with np.load(self._file_path, allow_pickle=False) as archive:
                self.data = archive['data']
        elif ext == '.bin':
            raw = np.fromfile(self._file_path, dtype=np.float64)
            if raw.size % 2:
                raise ValueError("binary single-scan data must contain m/z,intensity pairs")
            self.data = raw.reshape((-1, 2))
        else:
            raise ValueError(f"unsupported single-scan extension: {ext}")
        self.data = np.atleast_2d(self.data)
        if self.data.shape[1] < 2:
            raise ValueError("single-scan data must contain at least two columns")
        return self.data


    def get_all_scans(self):
        """Return a one-element list containing the m/z-intensity spectrum."""
        return [self.load_data()[:, :2]]

    def get_single_scan(self, scan=None):
        """Return the sole m/z-intensity spectrum; *scan* is ignored."""
        return self.load_data()[:, :2]

    def get_avg_scan(self, scan_range=None, time_range=None):
        """Return the sole spectrum; range arguments are ignored."""
        return self.load_data()[:, :2]

    def get_cdms_data(self):
        """Return five-column CD-MS events, filling absent metadata with defaults."""
        raw_dat = [self.load_data()]
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])

        try:
            scans = np.concatenate([d[:, 2] for d in raw_dat])
        except Exception as e:
            print("No scan data in NPZ file, populating with 1's")
            scans = np.ones_like(intensity)
        try:
            it = np.concatenate([d[:, 3] for d in raw_dat])
        except Exception as e:
            print("No injection time data in NPZ file, populating with 1's")
            it = np.ones_like(mz)

        try:
            times = np.concatenate([d[:, 4] for d in raw_dat])
        except Exception as e:
            print("No time data in NPZ file, populating with -1's")
            times = np.zeros_like(mz) - 1

        # elif ext.lower() == '.bin':
        #     try:
        #         raw_dat = raw_dat.reshape((int(len(raw_dat) / 3), 3))
        #     except Exception as e:
        #         raw_dat = raw_dat.reshape((int(len(raw_dat) / 2)), 2)
        #     mz = raw_dat[:, 0] > 0
        #     intensity = raw_dat[:, 1] > 0
        #     it = np.ones_like(mz)
        #     scans = np.ones_like(intensity)
        # else:
        #     it = np.ones_like(mz)
        #     scans = np.ones_like(intensity)

        return np.transpose([mz, intensity, scans, it, times])

    def get_imms_avg_scan(self, scan_range=None, time_range=None, mzbins=None):
        """Return the first three columns as m/z, drift time, and intensity."""
        del scan_range, time_range, mzbins
        self.immsdata = self.load_data()
        if self.immsdata.shape[1] < 3:
            raise ValueError("ion-mobility data must contain m/z, drift time, and intensity columns")
        self.immsdata = self.immsdata[:, :3]
        return self.immsdata

    def get_all_imms_scans(self):
        """Return a one-element list containing the ion-mobility array."""
        return [self.get_imms_avg_scan()]

    def get_imms_scan(self, s):
        """Return the sole ion-mobility scan; *s* is ignored."""
        return self.get_imms_avg_scan()


if __name__ == "__main__":
    path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\CDMS\\test_csv_cdms.csv"
    path = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Skippy\\HT SEC paper data\\Bgal GroEL\\20240412 Bgal GroEL bit5 zp10 3_2024-04-16-11-38-01_unidecfiles\\20240412 Bgal GroEL bit5 zp10 3_2024-04-16-11-38-01_rawdata.npz"
    path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\SingleScan\\test_csv.csv"
    path = r"Z:\Group Share\BHT\Q Exactive HF Data\RPLC-MS\Acquity UPLC\CDMS Injections to Stitch\20251216\20251216_BHT_0o1mgmL_carbonicanhydrase_CDMS_2uLinj_merged.npz"
    importer = SingleScanImporter(path)
    print(importer.get_cdms_data())







