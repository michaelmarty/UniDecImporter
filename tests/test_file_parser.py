from pathlib import Path

import h5py
import numpy as np
import pytest

from MassSpecImporter import FileParser


class FakeImporter:
    scans = np.array([1, 2, 3, 4])
    times = np.array([0.0, 1.0, 2.0, 3.0])

    def get_avg_scan(self, scan_range=None, time_range=None):
        del scan_range, time_range
        return np.array([[100.0, 2.0], [101.0, 3.0]])

    def get_max_time(self):
        return 3.0

    def get_max_scan(self):
        return 4

    def get_times_from_scans(self, scan_range):
        return [float(scan_range[0] - 1), np.mean(scan_range) - 1, float(scan_range[1] - 1)]


@pytest.fixture
def fake_source(tmp_path, monkeypatch):
    source = tmp_path / "sample.mzML"
    source.write_text("placeholder", encoding="utf-8")
    monkeypatch.setattr(FileParser.ImporterFactory, "create_importer", lambda path: FakeImporter())
    return source


def test_parse_text_and_missing_path(fake_source, tmp_path):
    outdir = tmp_path / "text-output"
    result = FileParser.parse(fake_source, [0, 1], 1, None, "slice", outdir)
    assert result == str(outdir / "slice_1.txt")
    np.testing.assert_array_equal(np.loadtxt(outdir / "slice_0.txt"), [[100, 2], [101, 3]])
    assert FileParser.parse(tmp_path / "missing", [0], 1, None, "x", outdir) is None
    with pytest.raises(ValueError):
        FileParser.parse(fake_source, [0], 1, None, "x", outdir, output="bad")


def test_parse_hdf5(fake_source, tmp_path):
    result = FileParser.parse(fake_source, [0, 1], 0.5, [10, 20], "slices", tmp_path, output="hdf5")
    with h5py.File(result) as hdf:
        assert hdf["ms_dataset"].attrs["num"] == 2
        assert hdf["ms_dataset/0"].attrs["Collision Voltage"] == 10
        np.testing.assert_array_equal(hdf["ms_dataset/1/raw_data"], [[100, 2], [101, 3]])


def test_parse_multiple(fake_source, tmp_path):
    second = tmp_path / "second.mzML"
    second.write_text("placeholder", encoding="utf-8")
    result = FileParser.parse_multiple([fake_source, second], 1, tmp_path, 0, 2,
                                       voltsarr=[[1, 2], [3, 4]], outputname="combined")
    with h5py.File(result) as hdf:
        assert hdf["ms_dataset"].attrs["num"] == 4
        assert hdf["ms_dataset/2"].attrs["Original File"] == "second.mzML"


def test_extract_time_and_scan_slices(fake_source, tmp_path):
    time_result = FileParser.extract(fake_source.name, tmp_path, timestep=1, output="txt")
    assert Path(time_result).exists()
    scan_result = FileParser.extract_scans(fake_source.name, tmp_path, scanbins=2, output="hdf5")
    with h5py.File(scan_result) as hdf:
        assert hdf["ms_dataset"].attrs["num"] == 2


def test_extract_timepoints_delegates(monkeypatch, tmp_path):
    seen = {}

    def fake_parse_multiple(*args):
        seen["args"] = args
        return "done"

    monkeypatch.setattr(FileParser, "parse_multiple", fake_parse_multiple)
    result = FileParser.extract_timepoints(["run_ramp_1_2_1.mzML"], [str(tmp_path)], 0, 2, 1, "out")
    assert result == "done"
    np.testing.assert_array_equal(seen["args"][5][0], [1, 2])


def test_extract_scans_multiple_files(fake_source, tmp_path):
    result = FileParser.extract_scans_multiple_files(
        [fake_source.name], [str(tmp_path)], startscan=1, endscan=3, outputname="multi"
    )
    with h5py.File(result) as hdf:
        assert hdf["ms_dataset"].attrs["num"] == 1
        assert hdf["ms_dataset/0"].attrs["Original File"] == fake_source.name


def test_get_files_is_case_insensitive(tmp_path, monkeypatch):
    (tmp_path / "one.mzML").write_text("", encoding="utf-8")
    (tmp_path / "two.MZML").write_text("", encoding="utf-8")
    (tmp_path / "ignore.txt").write_text("", encoding="utf-8")
    calls = []
    monkeypatch.setattr(FileParser, "extract", lambda *args, **kwargs: calls.append((args, kwargs)))
    FileParser.get_files(tmp_path, timestep=2, output="hdf5")
    assert len(calls) == 2
