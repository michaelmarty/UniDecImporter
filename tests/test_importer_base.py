import numpy as np
import pytest

from UniDecImporter.Importer import Importer


class DummyImporter(Importer):
    def __init__(self, path):
        super().__init__(path)
        self.scans = np.array([1, 3, 5])
        self.times = np.array([0.0, 1.0, 2.0])
        self.levels = np.array([1, 2, 1])
        self.scan_range = [1, 5]
        self.chrom_support = True
        self.cdms_support = True
        self.imms_support = True
        self._spectra = [
            np.array([[100, 1], [101, 2], [102, 1]], dtype=float),
            np.array([[100, 2], [101, 4], [102, 2]], dtype=float),
            np.array([[100, 3], [101, 6], [102, 3]], dtype=float),
        ]

    def get_all_scans(self):
        self.data = self._spectra
        return self.data

    def get_single_scan(self, scan):
        return self._spectra[self.get_scan_index(scan)]

    def get_avg_scan(self, scan_range=None, time_range=None):
        return self.avg_fast(scan_range, time_range)

    def get_all_imms_scans(self):
        self.immsdata = [
            np.array([[100, 1, 1], [100, 2, 2], [101, 1, 3], [101, 2, 4]], dtype=float)
            for _ in self.scans
        ]
        return self.immsdata


def test_base_metadata_and_lookup(tmp_path):
    path = tmp_path / "dummy.raw"
    path.write_bytes(b"x")
    importer = DummyImporter(path)
    assert importer.get_max_scan() == 5
    assert importer.get_max_time() == 2
    assert importer.get_scan_index(0) == 0
    assert importer.get_scan_index(4) == 2
    assert importer.get_scan_index(99) == 2
    assert importer.get_scan_time(3) == 1
    assert importer.get_time_scan(1.6) == 5
    assert importer.get_ms_order(3) == 2
    assert importer.get_scans_from_times((0.1, 1.9)) == [1, 5]
    assert importer.get_times_from_scans((1, 5)) == [0, 1, 2]
    with importer as entered:
        assert entered is importer


def test_base_ranges_averaging_and_centroiding(tmp_path):
    path = tmp_path / "dummy.raw"
    path.write_bytes(b"x")
    importer = DummyImporter(path)
    np.testing.assert_array_equal(importer.scan_range_from_inputs((-5, 99)), [1, 5])
    np.testing.assert_array_equal(importer.scan_range_from_inputs(time_range=(0.9, 0.9)), [3, 3])
    averaged = importer.get_avg_scan(scan_range=(1, 3))
    np.testing.assert_allclose(averaged[:, 1], [3, 6, 3])
    assert isinstance(importer.check_centroided(), bool)
    maxima = importer.get_mz_localmax(101, 10_000)
    assert maxima.shape == (3, 2)


def test_base_chromatogram_cdms_and_imms(tmp_path):
    path = tmp_path / "dummy.raw"
    path.write_bytes(b"x")
    importer = DummyImporter(path)
    importer.index_scans(0, 1)
    # Only MS1 scans are indexed.
    assert [scan.scan_num for scan in importer.indexed_file.indexed_scans] == [1, 5]
    assert importer.get_eic(101, 0.01).shape == (2, 2)
    assert importer.get_cdms_data().shape == (9, 5)
    assert importer.get_imms_scan(3).shape == (4, 3)
    assert importer.get_imms_avg_scan(mzbins=1).shape[1] == 3


def test_base_abstract_and_capability_errors(tmp_path):
    path = tmp_path / "base.txt"
    path.write_text("1 2\n", encoding="utf-8")
    base = Importer(path)
    with pytest.raises(NotImplementedError):
        base.get_all_scans()
    with pytest.raises(NotImplementedError):
        base.get_avg_scan()
    with pytest.raises(NotImplementedError):
        base.get_single_scan(1)
    with pytest.raises(Exception):
        base.get_tic()
    with pytest.raises(Exception):
        base.get_cdms_data()
    with pytest.raises(Exception):
        base.get_imms_avg_scan()
    with pytest.raises(FileNotFoundError):
        Importer(tmp_path / "missing")

