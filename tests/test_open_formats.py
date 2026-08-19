import re
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from UniDecImporter import get_importer
from UniDecImporter.MZML.mzML import (
    auto_gzip,
    correct_mzML_ID,
    gzip_files,
    get_data_from_spectrum as mzml_data,
    get_inj_time,
    search_by_id,
)
import UniDecImporter.MZML.mzML as mzml_module
from UniDecImporter.MZXML.mzXML import get_data_from_spectrum as mzxml_data


def require_real_data(path, minimum_size=1024):
    if not path.exists() or (path.is_file() and path.stat().st_size < minimum_size):
        pytest.skip(f"test data unavailable (run git lfs pull): {path}")
    return path


class FakeMzMLSpec:
    mz = np.array([5, 100, 100, 101], dtype=float)
    i = np.array([99, 2, 3, 0], dtype=float)


def test_mzml_spectrum_conversion_filtering_and_duplicates():
    data = mzml_data(FakeMzMLSpec(), check_duplicates=True)
    np.testing.assert_array_equal(data, [[100, 5], [101, 0]])
    np.testing.assert_array_equal(mzml_data(FakeMzMLSpec(), threshold=0), [[100, 2], [100, 3]])
    assert mzml_data(None) is None


def test_mzxml_spectrum_conversion():
    spectrum = {"m/z array": np.array([5, 100, 101]), "intensity array": np.array([9, 2, 0])}
    np.testing.assert_array_equal(mzxml_data(spectrum), [[100, 2], [101, 0]])
    np.testing.assert_array_equal(mzxml_data(spectrum, threshold=0), [[100, 2]])


def test_injection_time_parsing():
    root = ET.fromstring('<spectrum><cvParam name="ion injection time" value="12.5"/></spectrum>')
    spec = type("Spec", (), {"element": root})()
    assert get_inj_time(spec) == 12.5
    assert get_inj_time(type("Spec", (), {"element": ET.fromstring("<spectrum/>")})()) == 1


def test_search_by_id():
    obo = type("OBO", (), {
        "lookups": [{"MS:123": {"id": "MS:123", "name": "example"}}],
        "MS_tag_regex": re.compile(r"MS:\d+"),
    })()
    assert search_by_id(obo, 123) == "MS:123\nexample\n"
    assert search_by_id(obo, 999) == ""


def test_correct_mzml_id_preserves_input(tmp_path):
    source = tmp_path / "bad.mzML"
    original = '<mzML>\n<spectrum id="not-an-int" index="4">\n</spectrum>\n</mzML>\n'
    source.write_text(original, encoding="utf-8")
    corrected = correct_mzML_ID(source)
    assert source.read_text(encoding="utf-8") == original
    assert 'id="scan=5"' in open(corrected, encoding="utf-8").read()


def test_gzip_helpers_delegate_to_pymzml_indexer(tmp_path, monkeypatch):
    source = tmp_path / "tiny.mzML"
    source.write_text("<mzML/>", encoding="utf-8")
    calls = []
    fake_reader = type("Reader", (), {"get_spectrum_count": lambda self: 3})()
    monkeypatch.setattr(mzml_module.pymzml.run, "Reader", lambda path: fake_reader)
    monkeypatch.setattr(mzml_module, "index_gzip", lambda *args, **kwargs: calls.append((args, kwargs)))
    output = tmp_path / "tiny.mzML.gz"
    gzip_files(source, output)
    assert calls[0][0] == (source, output)
    assert calls[0][1]["max_idx"] == 35
    auto = auto_gzip(str(source))
    assert auto.endswith(".mzML.gz")
    assert len(calls) == 2


@pytest.fixture(scope="module", params=["test_mzml.mzML", "test_mzmlgz.mzML.gz", "test_mzxml.mzXML"])
def open_importer(request, data_dir):
    path = require_real_data(data_dir / request.param, minimum_size=1_000_000)
    importer = get_importer(path)
    yield importer
    importer.close()


@pytest.mark.integration
def test_open_format_metadata_and_single_scan(open_importer):
    importer = open_importer
    assert len(importer.scans) == len(importer.times) == len(importer.levels)
    assert len(importer.scans) > 20
    assert importer.get_max_scan() == importer.scans[-1]
    assert importer.get_max_time() == importer.times[-1]
    scan = importer.get_single_scan(importer.scans[0])
    assert scan.ndim == 2 and scan.shape[1] == 2 and len(scan) > 20
    assert importer.get_ms_order(importer.scans[0]) >= 1
    assert importer.get_polarity(importer.scans[0]) in {"Positive", "Negative"}


@pytest.mark.integration
def test_open_format_tic_average_and_eic(open_importer):
    importer = open_importer
    tic = importer.get_tic()
    assert tic.shape == (len(importer.scans), 2)
    average = importer.get_avg_scan(scan_range=(importer.scans[0], importer.scans[2]))
    assert average.ndim == 2 and average.shape[1] == 2 and len(average) > 20
    importer.index_scans(0, 1)
    eic = importer.get_eic(float(average[len(average) // 2, 0]), 0.5)
    assert eic.shape == (np.count_nonzero(importer.levels == 1), 2)


@pytest.mark.integration
def test_open_format_time_scan_roundtrip(open_importer):
    importer = open_importer
    start, end = importer.scans[1], importer.scans[-2]
    times = importer.get_times_from_scans((start, end))
    assert importer.get_scans_from_times((times[0], times[2])) == [start, end]


@pytest.mark.integration
def test_mzml_ion_mobility(data_dir):
    path = require_real_data(data_dir / "IMMS" / "test_agilentimms_mzml.mzML", 1_000_000)
    importer = get_importer(path)
    try:
        scans = importer.get_all_imms_scans()
        assert len(scans) > 1
        assert all(scan.shape[1] == 3 for scan in scans)
        one = importer.get_imms_scan(importer.scans[0])
        assert one.shape[1] == 3
        average = importer.get_imms_avg_scan(mzbins=1)
        assert average.shape[1] == 3
        assert importer.get_inj_time_array().shape == importer.scans.shape
        assert importer.get_property(importer.scans[0], "not present") == 1
        assert importer.get_isolation_mz_width(importer.scans[0]) == (None, None)
        cdms = importer.get_cdms_data()
        assert cdms.ndim == 2 and cdms.shape[1] == 5
    finally:
        importer.close()
