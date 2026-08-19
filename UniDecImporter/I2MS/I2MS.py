import sqlite3
import numpy as np
from ..Importer import Importer


class I2MSImporter(Importer):
    def __init__(self, file):
        super().__init__(file)
        self.conn = sqlite3.connect(file)
        self.cursor = self.conn.cursor()
        self.cursor.execute('SELECT * from main.Ion ORDER BY mz ASC')
        data = self.cursor.fetchall()
        self.data = np.array(data)
        self.cursor.execute('PRAGMA table_info(Ion)')
        keys = self.cursor.fetchall()
        self.keys = np.array(keys)
        self.mzkey = np.where(self.keys[:, 1] == "Mz")[0][0]
        self.scankey = np.where(self.keys[:, 1] == "ScanNumber")[0][0]
        self.slopekey = np.where(self.keys[:, 1] == "Slope")[0][0]
        self.scans = np.unique(self.data[:, self.scankey].astype(int))
        self.iitkey = np.where(self.keys[:, 1] == "InverseInjectionTimeSeconds")[0][0]
        self.invinjtime = self.data[:, self.iitkey]
        self.times = np.zeros(len(self.scans), dtype=float)
        self.scan_range = [int(np.min(self.scans)), int(np.max(self.scans))]
        '''
        try:
            self.invinjtime = np.where(self.keys[:, 1] == "InverseInjectionTimeSeconds")[0][0]
        except:
            self.invinjtime = None'''
        self.cdms_support=True
        self.imms_support=False
        self.chrom_support=False

    def get_all_scans(self, threshold=-1):
        scans = []
        for scan in self.scans:
            spectrum = self.get_single_scan(scan)
            if threshold >= 0:
                spectrum = spectrum[spectrum[:, 1] > threshold]
            scans.append(spectrum)
        return scans

    def get_cdms_data(self):
        mz = self.data[:, self.mzkey]
        intensity = self.data[:, self.slopekey]
        scans = self.data[:, self.scankey]
        if len(scans) != len(mz):
            scans = np.ones(len(mz))
        it = self.invinjtime
        times = np.zeros(len(mz)) - 1
        data_array = np.transpose([mz, intensity, scans, it, times])
        return data_array

    def get_cdms_data_by_scans(self, scan_range):
        data = self.get_cdms_data()
        mask = np.logical_and(data[:, 2] >= scan_range[0], data[:, 2] <= scan_range[1])
        return data[mask]

    def get_single_scan(self, scan=None):
        res = self.data[self.data[:, self.scankey] == scan]
        return np.transpose([res[:, self.mzkey], res[:, self.slopekey]])

    def close(self):
        if getattr(self, "cursor", None):
            self.cursor.close()
            self.cursor = None
        if getattr(self, "conn", None):
            self.conn.close()
            self.conn = None


    def get_scan_range(self):
        return list(self.scan_range)

    def get_avg_scan(self, bins=1, scan_range=None, time_range=None):
        del time_range
        if scan_range is not None:
            all_scans = self.get_cdms_data_by_scans(scan_range)[:, :2]
        else:
            all_scans = self.get_cdms_data()[:, :2]
        if len(all_scans) == 0:
            return np.empty((0, 2))
        width = float(bins)
        if width <= 0:
            raise ValueError("bins must be positive")
        mz, intensity = all_scans[:, 0], all_scans[:, 1]
        edges = np.arange(np.floor(mz.min() / width) * width,
                          np.ceil(mz.max() / width) * width + width, width)
        summed, edges = np.histogram(mz, bins=edges, weights=intensity)
        return np.column_stack((edges[:-1], summed))


if __name__ == "__main__":


    file = "Z:\\Group Share\\JGP\\DiverseDataExamples\\DataTypeCollection\\CDMS\\test_dmt_cdms.dmt"
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.use('WxAgg')
    i2ms = I2MSImporter(file)
    dat =i2ms.get_avg_scan()
    plt.plot(dat[:,0], dat[:,1])
    plt.show()
    exit()


    x = i2ms.data[:,7]
    y = i2ms.data[:,10]

    mzbins = 1000
    zbins = 0.1

    mzrange = [np.floor(np.amin(x)), np.amax(x)/5]
    zrange = [np.floor(np.amin(y)), np.amax(y)]
    mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + mzbins / 2, mzbins)
    # Weird fix to make this axis even is necessary for CuPy fft for some reason...
    if len(mzaxis) % 2 == 1:
        mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + 3 * mzbins / 2, mzbins)
    zaxis = np.arange(zrange[0] - zbins / 2., zrange[1] + zbins / 2, zbins)

    harray, xtab, ytab = np.histogram2d(x, y, [mzaxis, zaxis])
    xtab = xtab[1:] - mzbins / 2.
    ytab = ytab[1:] - zbins / 2.

    harray = np.transpose(harray)
    harray /= np.amax(harray)

    plt.imshow(harray, aspect="auto")
    plt.show()
