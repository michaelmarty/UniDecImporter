"""High-level Thermo RAW importer backed by RawFileReader."""

import numpy as np

from ..Importer import Importer
from .RawFileReader import RawFileReader as rr


class ThermoImporter(Importer):
    """
    Imports Thermo data files.

    Note: Thermo scans are 1 indexed, so the first scan is scan 1, not scan 0.
    """
    def __init__(self, path, silent=False, **kwargs):
        """Open a Thermo RAW file and initialize scan metadata."""
        super().__init__(path, **kwargs)

        if not silent:
            print("Launching Thermo Importer. If it fails, try to unblock the zip before unzipping")
        print("Reading Thermo Data:", path)
        self.msrun = rr(path)
        self.init_scans()

        self.cdms_support = True
        self.imms_support = False
        self.chrom_support = True
        self.thermo_support = True

    def init_scans(self):
        """Load the scan range and retention times from RawFileReader."""
        self.scan_range = self.msrun.scan_range
        self.scans = np.arange(self.scan_range[0], self.scan_range[1] + 1)
        self.times = []
        for s in self.scans:
            self.times.append(self.msrun.scan_time_from_scan_name(s))

        self.times = np.array(self.times)
        print("Number of Scans", len(self.scans))

    def get_all_scans(self, threshold=-1):
        """Load all profile scans, optionally filtering by intensity."""
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.GetSpectrum(s))
            impdat = impdat[impdat[:, 0] > 10]
            if threshold >= 0:
                impdat = impdat[impdat[:, 1] > threshold]
            self.data.append(impdat)
        return self.data

    def get_all_centroid_scans(self, threshold=-1):
        """Load all centroid scans, optionally filtering by intensity."""
        self.data = []
        for s in self.scans:
            impdat = np.array(self.msrun.GetCentroidSpectrum(s))
            impdat = impdat[impdat[:, 0] > 10]
            if threshold >= 0:
                impdat = impdat[impdat[:, 1] > threshold]
            self.data.append(impdat)
        return self.data

    def get_single_scan(self, s):
        """Return one profile scan as an m/z-intensity array."""
        impdat = np.array(self.msrun.GetSpectrum(int(s)))
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def grab_centroid_data(self, s):
        """Return one centroid scan as an m/z-intensity array."""
        impdat = np.array(self.msrun.GetCentroidSpectrum(int(s)))
        impdat = impdat[impdat[:, 0] > 10]
        return impdat

    def grab_all_centroid_dat(self):
        """Concatenate centroid peaks from every scan."""
        data = []
        for i in range(len(self.scans)):
            impdat = np.array(self.msrun.GetCentroidSpectrum(self.scans[i]))

            impdat = impdat[impdat[:, 0] > 10]
            data.append(impdat)
        if len(data) > 0:
            return np.vstack(data)
        else:
            return np.empty((0, 2))

    def get_avg_scan(self, scan_range=None, time_range=None):
        """Return a RawFileReader average over an inclusive scan or time range."""
        scan_range = self.scan_range_from_inputs(scan_range, time_range)

        if scan_range[1] - scan_range[0] > 1:
            print("Getting Data from Scans:", scan_range)
            scan_range = [scan_range[0], scan_range[1]]
            data = np.array(list(self.msrun.GetAverageSpectrum(scan_range)))
        else:
            print("Getting Data from Scan:", scan_range[0])
            impdat = self.get_single_scan(scan_range[0])
            data = impdat

        return data

    def get_tic(self):
        """Return the total-ion chromatogram."""
        return self.msrun.GetChromatogram()

    def get_eic(self, mass, mz_tol, rt_range=None):
        """Return an extracted-ion chromatogram for the selected time range."""
        mass_range = [mass-(mz_tol/2), mass+(mz_tol/2)]
        if rt_range is not None:
            min_idx = int(np.argmin(np.abs(self.times - rt_range[0])))
            max_idx = int(np.argmin(np.abs(self.times - rt_range[1])))
            scan_range = [self.scans[min_idx], self.scans[max_idx]]
        else:
            scan_range = [self.scans[0], self.scans[-1]]
        return self.msrun.Get_EIC(massrange=mass_range, scanrange=scan_range)

    def get_inj_time_array(self):
        """Return ion injection times aligned with the scan list."""
        its = []
        for i, s in enumerate(self.scans):
            it, res, an1, an2 = self.msrun.get_scan_header(s)
            try:
                it = float(it)
            except:
                print("Error in scan header:", i, s, it)
                it = 1
            its.append(it)
        return np.array(its)

    def get_analog_voltage1(self):
        """Return the first analog-voltage header value for each scan."""
        vs = []
        for i, s in enumerate(self.scans):
            it, res, an1, an2 = self.msrun.get_scan_header(s)
            vs.append(an1)
        return np.array(vs)

    def get_analog_voltage2(self):
        """Return the second analog-voltage header value for each scan."""
        vs = []
        for i, s in enumerate(self.scans):
            it, res, an1, an2 = self.msrun.get_scan_header(s)
            vs.append(an2)
        return np.array(vs)

    def get_polarity(self, scan=1):
        """Return the polarity encoded in a scan event string."""
        # print(dir(self.msrun.source))
        scan_mode = self.msrun.source.GetScanEventStringForScanNumber(int(scan))
        if "+" in scan_mode:
            print("Polarity: Positive")
            return "Positive"
        if "-" in scan_mode[:10]:
            print("Polarity: Negative")
            return "Negative"
        print("Polarity: Unknown")
        return None

    def get_ms_order(self, scan=1):
        """Return the MS level reported by RawFileReader."""
        order = self.msrun.GetMSOrder(int(scan))
        return order

    def get_isolation_mz_width(self, s):
        """Return precursor m/z and isolation width for a reaction scan."""
        scanFilter = self.msrun.GetScanFilter(int(s))
        reaction = scanFilter.GetReaction(0)
        mz = reaction.PrecursorMass
        width = reaction.IsolationWidth
        return mz, width

    def get_cdms_data(self, scan_range=None):
        """Return centroid peaks as injection-corrected five-column CD-MS events."""
        raw_dat = self.get_all_centroid_scans(threshold=0)
        scans = self.scans

        it = 1. / self.get_inj_time_array()
        mz = np.concatenate([d[:, 0] for d in raw_dat])
        scans = np.concatenate([s * np.ones(len(raw_dat[i])) for i, s in enumerate(self.scans)])
        try:
            intensity = np.concatenate([d[:, 1] * it[i] / 1000. for i, d in enumerate(raw_dat)])
        except Exception as e:
            print("Mark1:", e, it)
            intensity = np.concatenate([d[:, 1] for i, d in enumerate(raw_dat)])
        try:
            it = np.concatenate([it * np.ones(len(raw_dat[i])) for i, it in enumerate(it)])
        except Exception as e:
            print("Error with injection time correction:", e)

        try:
            times = np.concatenate([self.times[i] * np.ones(len(raw_dat[i])) for i, _ in enumerate(raw_dat)])
        except Exception as e:
            print("Error with time array:", e)
            times = np.zeros_like(mz) - 1

        data_array = np.transpose([mz, intensity, scans, it, times])
        return data_array

    def close(self):
        """Close the Thermo RawFileReader handle."""
        self.msrun.Close()
        return


if __name__ == "__main__":
    # import matplotlib.pyplot as plt
    test = "C:\\Python\\UniDec3\\TestSpectra\\test.raw"
    test = "Z:\\Group Share\\Annika\\2025-04-29_OE240_Stellar_UTAustin\\Analytical Flow\\Stellar\\Data\\2025-03-13_Stellar_UTAustin_Sample_E_posneg_tMS2_01.raw"
    d = ThermoImporter(test, silent=False)
    data = d.get_single_scan(100)
    print(data)
    # plt.plot(data[:, 0], data[:, 1])
    # plt.show()
    exit()
    cdms_dat = importer.get_cdms_data()
    for i in cdms_dat:
        print(i)
