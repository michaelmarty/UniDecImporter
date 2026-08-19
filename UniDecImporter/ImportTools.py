import re
import numpy as np
from collections import defaultdict
from . import tools as ud
from copy import deepcopy
import math

def header_test(path, deletechars=None, delimiter=" |\t|,", strip_end_space=True):
    """
    A quick test of a file to remove any problems from the header.

    If scans through each line of an input file specificed by path and count the number of lines that cannot be
    completely converted to floats. It breaks the loop on the first hit that is entirely floats, so any header lines
    blow that will not be counted.
    :param path: Path of file to be scanned
    :param deletechars: Characters to be removed from the file before scanning
    :param delimiter: Delimiter to be used in the file
    :param strip_end_space: Boolean to strip the end space from the line
    :return: Integer length of the number of header lines
    """
    header = 0
    try:
        with open(path, "r", encoding="utf-8-sig") as f:
            for line in f:
                # If the line is empty, skip it and add to the header
                if line == "\n":
                    header += 1
                    continue
                # Remove any characters that are defined in deletechars
                if deletechars is not None:
                    for c in deletechars:
                        line = line.replace(c, "")

                # Strip the end space if strip_end_space is True
                if strip_end_space:
                    if line[-1] == "\n" and len(line) > 1:
                        line = line[:-1]
                    if line[-1] == " " and len(line) > 1:
                        line = line[:-1]

                numeric = True
                # Split the line by the delimiter and try to convert each element to a float
                for sin in re.split(delimiter, line):
                    if sin == "":
                        continue
                    try:
                        float(sin)
                    except ValueError:
                        header += 1
                        numeric = False
                        break
                if numeric:
                    break
        # If the header is greater than 0, print the header length
        if header > 0:
            print("Header Length:", header)
    except (ImportError, OSError, AttributeError, IOError) as e:
        print("Failed header test", e)
        header = 0
    return int(header)


def get_resolution(testdata):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    testdata = np.asarray(testdata)
    if testdata.ndim != 2 or testdata.shape[1] < 2 or len(testdata) < 2:
        raise ValueError("at least two m/z-intensity rows are required")
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    resolutions = ud.safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)


def fit_line(x, a, b):
    return a * x ** b


def get_longest_index(datalist):
    if len(datalist) == 0:
        raise ValueError("datalist must contain at least one spectrum")
    lengths = [len(x) for x in datalist]
    return np.argmax(lengths)

def merge_spectra(datalist, mzbins=None, type="Integrate"):
    """
    Merge together a list of data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Filter out junk spectra that are empty
    datalist = [x for x in datalist if len(x) > 0]
    if not datalist:
        return np.empty((0, 2))
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    concat = np.concatenate(datalist)
    # xvals = concat[:, 0]
    # print "Median Resolution:", resolution
    # axis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    # for d in datalist:
    #    print(d)
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.

    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution(datalist[maxlenpos])
        if resolution < 0:
            print("ERROR with auto resolution:", resolution, maxlenpos, datalist[maxlenpos])
            print("Using ABS")
            resolution = np.abs(resolution)
        elif resolution == 0:
            print("ERROR, resolution is 0, using 20000.", maxlenpos, datalist[maxlenpos])
            resolution = 20000

        axis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        if float(mzbins) < 0:
            raise ValueError("mzbins must be non-negative")
        axis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]) + float(mzbins) / 2,
                         float(mzbins))
    template = np.transpose([axis, np.zeros_like(axis)])

    print("Length merge axis:", len(template))

    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = ud.mergedata(template, d)
            elif type == "Integrate":
                newdat = ud.lintegrate(d, axis, fastmode=True)
            else:
                print("ERROR: unrecognized merge spectra type:", type)
                continue

            template[:, 1] += newdat[:, 1]

    # Trying to catch if the data is backwards
    try:
        if template[1, 0] < template[0, 0]:
            template = template[::-1]
    except:
        pass
    return template

def get_resolution_im(data):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    testdata = deepcopy(data)
    testdata = testdata[testdata[:, 2] > 0]
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    b1 = diffs[:,1] > 0
    diffs = diffs[b1]
    resolutions = ud.safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)

def merge_im_spectra(datalist, mzbins=None, type="Integrate"):
    """
    Merge together a list of ion mobility data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    datalist = [np.asarray(x) for x in datalist if len(x) > 0]
    if not datalist:
        return np.empty((0, 3))
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    if len(datalist) > 1:
        concat = np.concatenate(datalist)
    else:
        concat = np.array(datalist[0])
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.
    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution_im(datalist[maxlenpos])
        mzaxis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        if float(mzbins) < 0:
            raise ValueError("mzbins must be non-negative")
        mzaxis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]) + float(mzbins) / 2,
                           float(mzbins))

    # For drift time, use just unique drift time values. May need to make this fancier.
    dtaxis = np.sort(np.unique(concat[:, 1]))

    # Create the mesh grid from the new axes
    X, Y = np.meshgrid(mzaxis, dtaxis, indexing="ij")

    template = np.transpose([np.ravel(X), np.ravel(Y), np.ravel(np.zeros_like(X))])
    print("Shape merge axis:", X.shape)
    if len(mzaxis) == 1:
        half_width = float(mzbins or 1) / 2
        xbins = np.array([mzaxis[0] - half_width, mzaxis[0] + half_width])
    else:
        xdiff = np.diff(mzaxis)
        xbins = np.concatenate(([mzaxis[0] - xdiff[0] / 2],
                                mzaxis[:-1] + xdiff / 2,
                                [mzaxis[-1] + xdiff[-1] / 2]))
    if len(dtaxis) == 1:
        ybins = np.array([dtaxis[0] - 0.5, dtaxis[0] + 0.5])
    else:
        ydiff = np.diff(dtaxis)
        ybins = np.concatenate(([dtaxis[0] - ydiff[0] / 2],
                                dtaxis[:-1] + ydiff / 2,
                                [dtaxis[-1] + ydiff[-1] / 2]))
    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                raise ValueError("2D interpolation is not supported; use type='Integrate'")
            elif type == "Integrate":
                newdat, xedges, yedges = np.histogram2d(d[:, 0], d[:, 1], bins=[xbins, ybins], weights=d[:, 2])
            else:
                print("ERROR: unrecognized merge spectra type:", type)
                continue
            template[:, 2] += np.ravel(newdat)
    return template


def nonlinear_axis(start, end, res):
    """
    Creates a nonlinear axis with the m/z values spaced with a defined and constant resolution.
    :param start: Minimum m/z value
    :param end: Maximum m/z value
    :param res: Resolution of the axis ( m / delta m)
    :return: One dimensional array of the nonlinear axis.
    """
    axis = []
    i = start
    axis.append(i)
    i += i / fit_line(i, res[0], res[1])
    while i < end:
        axis.append(i)
        i += i / fit_line(i, res[0], res[1])
    return np.array(axis)

#
# def waters_convert2(path, config=None, outfile=None, time_range=None):
#     data = WDI(path).get_data(time_range=time_range)
#
#     if outfile is None:
#         outfile = os.path.join(path, "converted_rawdata.txt")
#     np.savetxt(outfile, data)
#     return data

def compute_bin_indices(mzs, intensities, min_mz, bin_width):
    if bin_width <= 0:
        raise ValueError("bin_width must be positive")
    indices = ((mzs - min_mz) / bin_width).astype(np.int32)
    return indices, mzs, intensities

class IndexedScan:
    """
    A class to store an in individual scan with the peaks indexed into bins by m/z.
    """
    def __init__(self, scan, rt, scan_num, min_mz, bin_width):
        self.scan = scan
        self.retention_time = rt
        self.scan_num = scan_num
        self.min_mz = min_mz
        self.bin_width = bin_width
        self.indexed_peaks = self.index_peaks(self.min_mz, self.bin_width)


    def index_peaks(self, min_mz, bin_width):
        mzs = np.array(self.scan[:,0])
        intensities = np.array(self.scan[:,1])
        indices, mzs, intensities = compute_bin_indices(mzs, intensities, min_mz, bin_width)

        indexed_peaks = defaultdict(list)
        for mz, intensity, idx in zip(mzs, intensities, indices):
            indexed_peaks[idx].append([mz, intensity])

        return indexed_peaks

    def get_intensity(self, mz, mz_tolerance):
        """
        Calculates the total intensity within :param mz_tolerance
        around :param mz within the scan.
        """
        #this is easy if the mz+/- the tolerance is within a single bin..
        min_index = math.trunc((mz - mz_tolerance - self.min_mz) / self.bin_width)
        max_index = math.trunc((mz + mz_tolerance - self.min_mz) / self.bin_width)
        if min_index == max_index:
            sum_intensity = 0
            #sum the intensities of all peaks in the bin within the tolerance of the mz
            if self.indexed_peaks.get(min_index, None) is None:
                return 0
            else:
                bin_peaks = self.indexed_peaks.get(min_index)
                for i in range(len(bin_peaks)):
                    if abs(bin_peaks[i][0] - mz) <= mz_tolerance:
                        sum_intensity += bin_peaks[i][1]
            return sum_intensity
        else:
            #check peaks from both bins..
            sum_intensity = 0
            indices = range(min_index, max_index + 1, 1)
            for i in range(len(indices)):
                if self.indexed_peaks.get(indices[i], None) is None:
                    continue
                else:
                    bin_peaks = self.indexed_peaks.get(indices[i])
                    for j in range(len(bin_peaks)):
                        if abs(bin_peaks[j][0] - mz) <= mz_tolerance:
                            sum_intensity += bin_peaks[j][1]
            return sum_intensity

    def get_exp_mz(self, mz, mz_tol):
        """
        Takes an m/z and an absolute width and finds the nearest experimental m/z in the scan
        within the given tolerance to the query m/z.

        :param mz: the target m/z
        :param mz_tol: the width (in m/z space) within which
        """
        min_index = math.trunc((mz - mz_tol - self.min_mz) / self.bin_width)
        max_index = math.trunc((mz + mz_tol - self.min_mz) / self.bin_width)
        if min_index == max_index:
            peaks_sorted = sorted(self.indexed_peaks.get(min_index, []), key=lambda x:x[1], reverse=True)
        else:
            indices = range(min_index, max_index+1, 1)
            peaks_within_tol = []
            for s in indices:
                if self.indexed_peaks.get(s, None) is None:
                    continue
                else:
                    bin_peaks = self.indexed_peaks.get(s)
                    for p in bin_peaks:
                        if abs(p[0] - mz) <= mz_tol:
                            peaks_within_tol.append(p)
            peaks_sorted = sorted(peaks_within_tol, key=lambda x:x[1], reverse=True)

        if not peaks_sorted:
            return None
        else:
            return peaks_sorted[0][0]

class IndexedFile:
    """
    A class to store a set of indexed scans originating from a single mzML file.
    """
    def __init__(self):
        self.min_mz = 0
        self.bin_width = 1
        self.indexed_scans = []

    def extract_xic(self, mz, mz_width, rt_range=None):
        """
        Takes an mz and an mz width (and optionally an rt range)
        and produces an eic (scans, retention times, and intensities)
        """
        rts = []
        scans = []
        intensities = []

        t_min = -1
        t_max = -1
        if rt_range is not None:
            t_min = rt_range[0]
            t_max = rt_range[1]

        for i in range(len(self.indexed_scans)):
            if rt_range is not None and (self.indexed_scans[i].retention_time < t_min or self.indexed_scans[i].retention_time > t_max):
                continue
            elif rt_range is None or (t_min < self.indexed_scans[i].retention_time < t_max):
                intensity = self.indexed_scans[i].get_intensity(mz, mz_width)
                scans.append(self.indexed_scans[i].scan_num)
                rts.append(self.indexed_scans[i].retention_time)
                intensities.append(intensity)
        return np.transpose([rts, intensities])

    def get_indexed_spectrum_atRt(self, rt):
        """
        Grabs the spectrum with an RT nearest to the provided rt (minutes)
        """
        if not self.indexed_scans:
            raise LookupError("no scans have been indexed")
        index = int(np.argmin([abs(scan.retention_time - rt) for scan in self.indexed_scans]))
        return self.indexed_scans[index]
