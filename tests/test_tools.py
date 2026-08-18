import numpy as np
import pytest

from MassSpecImporter import tools
from MassSpecImporter.ImportTools import (
    IndexedFile,
    IndexedScan,
    compute_bin_indices,
    fit_line,
    get_longest_index,
    get_resolution,
    get_resolution_im,
    header_test,
    merge_im_spectra,
    merge_spectra,
    nonlinear_axis,
)


def test_core_numerical_helpers():
    assert tools.isempty(None)
    assert tools.isempty([])
    assert not tools.isempty([1])
    np.testing.assert_allclose(tools.safedivide(np.array([4, 3]), np.array([2, 0])), [2, 0])
    assert tools.nearestunsorted([10, 4, 8], 7.5) == 2
    np.testing.assert_allclose(tools.nonlinear_axis(100, 103, 100), [100, 101, 102.01])
    with pytest.raises(ValueError):
        tools.nonlinear_axis(0, 10, 100)


def test_interpolation_binning_and_peak_extraction():
    template = np.column_stack(([1, 2, 3], np.zeros(3)))
    data = np.array([[1, 2], [3, 6]])
    np.testing.assert_allclose(tools.mergedata(template, data)[:, 1], [2, 4, 6])

    binned = tools.lintegrate(np.array([[0.9, 2], [1.1, 3], [2.1, 5]]), np.array([1, 2]))
    assert binned[:, 1].sum() == 10
    assert tools.data_extract(np.array([[9.9, 2], [10.0, 7], [10.1, 3]]), 10, 1, 0.2) == 7
    assert tools.data_extract(np.array([[9.9, 2], [10.0, 7]]), 10, 4, 0.2) == 10
    assert tools.data_extract(np.empty((0, 2)), 10, 1, 1) == 0
    with pytest.raises(ValueError):
        tools.data_extract(data, 2, 99)


def test_autocorrelation_ratio():
    assert tools.get_autocorr_ratio(np.array([[1, 0], [2, 0]])) == 0
    assert tools.get_autocorr_ratio(np.array([[1, 1], [2, 1], [3, 1]])) == pytest.approx(2 / 3)


def test_header_test_stops_at_first_numeric_line(tmp_path):
    path = tmp_path / "header.csv"
    path.write_text("title\nunits,counts\n1,2\nnot,data\n", encoding="utf-8")
    assert header_test(path) == 2
    assert header_test(tmp_path / "missing.txt") == 0


def test_resolution_and_axis_helpers():
    spectrum = np.array([[100, 1], [101, 2], [102, 3]], dtype=float)
    assert get_resolution(spectrum) == pytest.approx(101.5)
    with pytest.raises(ValueError):
        get_resolution(spectrum[:1])
    assert fit_line(4, 2, 0.5) == 4
    assert get_longest_index([np.zeros((2, 2)), np.zeros((3, 2))]) == 1
    with pytest.raises(ValueError):
        get_longest_index([])
    assert len(nonlinear_axis(100, 102, (100, 0))) == 2


def test_merge_spectra_integrates_and_interpolates():
    first = np.array([[100, 1], [101, 2], [102, 3]], dtype=float)
    second = np.array([[100, 3], [101, 2], [102, 1]], dtype=float)
    integrated = merge_spectra([first, second, np.empty((0, 2))], mzbins=1)
    np.testing.assert_allclose(integrated[:, 1], [4, 4, 4])
    interpolated = merge_spectra([first, second], mzbins=1, type="Interpolate")
    np.testing.assert_allclose(interpolated[:, 1], [4, 4, 4])
    assert merge_spectra([], mzbins=1).shape == (0, 2)
    with pytest.raises(ValueError):
        merge_spectra([first], mzbins=-1)


def test_merge_ion_mobility_spectra():
    data = np.array([
        [100, 1, 2], [100, 2, 3], [101, 1, 4], [101, 2, 5],
    ], dtype=float)
    assert get_resolution_im(data) > 0
    merged = merge_im_spectra([data, data], mzbins=1)
    assert merged.shape == (4, 3)
    assert merged[:, 2].sum() == 2 * data[:, 2].sum()
    assert merge_im_spectra([], mzbins=1).shape == (0, 3)
    with pytest.raises(ValueError):
        merge_im_spectra([data], mzbins=1, type="Interpolate")


def test_indexed_scan_and_file():
    spectrum = np.array([[99.9, 2], [100.05, 5], [101.0, 3]])
    indices, mzs, intensities = compute_bin_indices(spectrum[:, 0], spectrum[:, 1], 0, 1)
    assert list(indices) == [99, 100, 101]
    np.testing.assert_array_equal(mzs, spectrum[:, 0])
    np.testing.assert_array_equal(intensities, spectrum[:, 1])
    with pytest.raises(ValueError):
        compute_bin_indices(mzs, intensities, 0, 0)

    scan = IndexedScan(spectrum, rt=2.0, scan_num=7, min_mz=0, bin_width=1)
    assert scan.get_intensity(100, 0.11) == 7
    assert scan.get_intensity(500, 0.1) == 0
    assert scan.get_exp_mz(100, 0.2) == 100.05
    assert scan.get_exp_mz(500, 0.1) is None

    indexed = IndexedFile()
    indexed.indexed_scans = [scan, IndexedScan(spectrum, 4.0, 8, 0, 1)]
    xic = indexed.extract_xic(100, 0.11)
    np.testing.assert_allclose(xic, [[2, 7], [4, 7]])
    assert indexed.extract_xic(100, 0.11, (2.5, 4.5)).shape == (1, 2)
    assert indexed.get_indexed_spectrum_atRt(3.8).scan_num == 8
    with pytest.raises(LookupError):
        IndexedFile().get_indexed_spectrum_atRt(1)

