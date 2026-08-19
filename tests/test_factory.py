import importlib
import platform

import pytest

import UniDecImporter
factory_module = importlib.import_module("UniDecImporter.ImporterFactory")
from UniDecImporter import ImporterFactory, get_importer
from UniDecImporter.I2MS.I2MS import I2MSImporter
from UniDecImporter.MZML.mzML import MZMLImporter
from UniDecImporter.MZXML.mzXML import MZXMLImporter
from UniDecImporter.SingleScanImporter import SingleScanImporter
from UniDecImporter.exceptions import UnsupportedFormatError, VendorReaderUnavailableError


def test_public_api_and_factory_instance():
    assert UniDecImporter.__version__ == "0.1.0"
    assert get_importer is not None
    assert ".mzml" in ImporterFactory().recognized_file_types
    assert ".raw" in UniDecImporter.recognized_types


@pytest.mark.integration
@pytest.mark.parametrize(
    ("relative", "expected"),
    [
        ("SingleScan/test_csv.csv", SingleScanImporter),
        ("CDMS/test_dmt_cdms.dmt", I2MSImporter),
        ("test_mzml.mzML", MZMLImporter),
        ("test_mzmlgz.mzML.gz", MZMLImporter),
        ("test_mzxml.mzXML", MZXMLImporter),
    ],
)
def test_factory_routes_open_formats(data_dir, relative, expected):
    importer = get_importer(data_dir / relative)
    try:
        assert isinstance(importer, expected)
    finally:
        importer.close()


def test_factory_rejects_unknown_extension(tmp_path):
    path = tmp_path / "sample.nope"
    path.write_text("data", encoding="utf-8")
    with pytest.raises(UnsupportedFormatError, match=".nope"):
        get_importer(path)


@pytest.mark.parametrize(("name", "vendor"), [("sample.RAW", "Thermo"), ("sample.d", "Agilent")])
def test_vendor_format_has_explicit_platform_error(tmp_path, monkeypatch, name, vendor):
    path = tmp_path / name
    path.write_text("data", encoding="utf-8")
    monkeypatch.setattr(factory_module, "_is_windows_x64", lambda: False)
    with pytest.raises(VendorReaderUnavailableError, match=vendor):
        get_importer(path)


def test_waters_directory_has_explicit_platform_error(tmp_path, monkeypatch):
    path = tmp_path / "sample.raw"
    path.mkdir()
    monkeypatch.setattr(factory_module, "_is_windows_x64", lambda: False)
    with pytest.raises(VendorReaderUnavailableError, match="Waters"):
        get_importer(path)


def test_agilent_reports_missing_runtime_on_windows_x64(tmp_path, monkeypatch):
    path = tmp_path / "sample.d"
    path.mkdir()
    monkeypatch.setattr(factory_module, "_is_windows_x64", lambda: True)
    monkeypatch.delenv("AGILENT_DA_SDK_DIR", raising=False)
    with pytest.raises(VendorReaderUnavailableError, match="Agilent reader"):
        get_importer(path)


def test_platform_predicate_matches_runtime():
    expected = platform.system() == "Windows" and platform.machine().lower() in {"amd64", "x86_64", "x64"}
    assert factory_module._is_windows_x64() is expected
