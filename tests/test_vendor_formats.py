import platform

import numpy as np
import pytest

from UniDecImporter import get_importer
from UniDecImporter.ImporterFactory import VendorReaderUnavailableError


WINDOWS_X64 = platform.system() == "Windows" and platform.machine().lower() in {"amd64", "x86_64", "x64"}


def open_vendor_or_skip(path, **kwargs):
    if not WINDOWS_X64:
        pytest.skip("vendor readers require Windows x86-64")
    try:
        return get_importer(path, **kwargs)
    except VendorReaderUnavailableError as error:
        pytest.skip(str(error))


@pytest.mark.vendor
@pytest.mark.integration
def test_thermo_lcms_reader(data_dir):
    importer = open_vendor_or_skip(data_dir / "test_thermo.RAW", silent=True)
    try:
        assert importer.thermo_support and importer.chrom_support and importer.cdms_support
        assert len(importer.scans) > 20
        scan = importer.get_single_scan(importer.scans[0])
        assert scan.ndim == 2 and scan.shape[1] == 2
        average = importer.get_avg_scan(scan_range=(importer.scans[0], importer.scans[2]))
        assert average.shape[1] == 2
        assert importer.get_tic().shape[1] == 2
        assert importer.get_eic(1000, 1).shape[1] == 2
        assert importer.get_inj_time_array().shape == importer.scans.shape
        assert importer.get_analog_voltage1().shape == importer.scans.shape
        assert importer.get_analog_voltage2().shape == importer.scans.shape
        assert importer.get_polarity(importer.scans[0]) in {"Positive", "Negative"}
        assert importer.get_ms_order(importer.scans[0]) >= 1
    finally:
        importer.close()


@pytest.mark.vendor
@pytest.mark.integration
def test_thermo_cdms_reader(data_dir):
    importer = open_vendor_or_skip(data_dir / "CDMS" / "test_raw_cdms.RAW", silent=True)
    try:
        centroid = importer.grab_centroid_data(importer.scans[0])
        assert centroid.shape[1] == 2
        assert importer.grab_all_centroid_dat().shape[1] == 2
        cdms = importer.get_cdms_data()
        assert cdms.ndim == 2 and cdms.shape[1] == 5 and len(cdms) > 100
    finally:
        importer.close()


@pytest.mark.vendor
@pytest.mark.integration
def test_waters_lcms_reader(data_dir):
    importer = open_vendor_or_skip(data_dir / "test_waters.raw")
    try:
        assert importer.chrom_support and importer.imms_support
        assert len(importer.scans) > 20
        assert importer.get_single_scan(1).shape[1] == 2
        scans = importer.get_all_scans()
        assert len(scans) == len(importer.scans)
        assert importer.get_avg_scan(scan_range=(1, 3)).shape[1] == 2
        assert importer.get_tic().shape == (len(importer.scans), 2)
        assert importer.get_bpi().shape == (len(importer.scans), 2)
        assert importer.get_polarity() in {"Positive", "Negative"}
        assert importer.get_ms_order() >= 1
        importer.index_scans(0, 1)
        assert importer.get_eic(1000, 1).shape == (len(importer.scans), 2)
    finally:
        importer.close()


@pytest.mark.vendor
@pytest.mark.integration
def test_waters_ion_mobility_reader(data_dir):
    importer = open_vendor_or_skip(data_dir / "IMMS" / "test_waters.raw")
    try:
        scan = importer.get_imms_scan(1)
        assert scan.ndim == 2 and scan.shape[1] == 3 and len(scan) > 20
        all_scans = importer.get_all_imms_scans()
        assert len(all_scans) == len(importer.scans)
        average = importer.get_imms_avg_scan(scan_range=(1, 2), mzbins=1)
        assert average.shape[1] == 3
        importer.get_stats()
        assert importer.get_stat_code(int(importer.stat_nums[0])) is not None
        assert importer.get_stat_name(str(importer.stat_names[0])) is not None
    finally:
        importer.close()


@pytest.mark.vendor
@pytest.mark.integration
def test_agilent_lcms_reader(data_dir):
    importer = open_vendor_or_skip(data_dir / "test_agilent.d")
    try:
        assert importer.chrom_support and not importer.cdms_support and not importer.imms_support
        assert len(importer.scans) > 1
        scan = importer.get_single_scan(importer.scans[0])
        assert scan.ndim == 2 and scan.shape[1] == 2
        assert importer.get_avg_scan(scan_range=(importer.scans[0], importer.scans[2])).shape[1] == 2
        assert importer.get_tic().shape[1] == 2
        assert importer.get_eic(1000, 1).shape[1] == 2
        assert importer.get_polarity(importer.scans[0]) in {"Positive", "Negative"}
        assert importer.get_ms_order(importer.scans[0]) >= 1
    finally:
        importer.close()


def test_vendor_environment_variables_are_documented_names():
    # A regression check for the deployment hooks; no SDK is loaded here.
    assert "THERMO_RAW_FILE_READER_DIR" not in {"", None}
    assert "MASSLYNX_RAW_DLL" not in {"", None}
    assert "AGILENT_DA_SDK_DIR" not in {"", None}
