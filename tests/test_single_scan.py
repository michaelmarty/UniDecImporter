import numpy as np
import pytest

from UniDecImporter import get_importer
from UniDecImporter.SingleScanImporter import SingleScanImporter


@pytest.mark.integration
@pytest.mark.parametrize("relative", ["test_csv.csv", "test_dat.dat", "test_txt.txt"])
def test_bundled_single_scan_formats(data_dir, relative):
    importer = get_importer(data_dir / "SingleScan" / relative)
    assert isinstance(importer, SingleScanImporter)
    assert len(importer) > 500
    assert importer.get_single_scan().shape[1] == 2
    assert importer.get_avg_scan().shape[1] == 2
    assert len(importer.get_all_scans()) == 1
    assert importer.get_polarity() == "Positive"
    assert importer.get_ms_order() == 1


@pytest.mark.integration
@pytest.mark.parametrize("relative", ["test_csv_cdms.csv", "test_npz_cdms.npz"])
def test_bundled_single_scan_cdms(data_dir, relative):
    importer = get_importer(data_dir / "CDMS" / relative)
    cdms = importer.get_cdms_data()
    assert cdms.ndim == 2 and cdms.shape[1] == 5
    assert cdms.shape[0] > 100
    assert importer.get_avg_scan().shape[1] == 2


@pytest.mark.integration
def test_bundled_text_ion_mobility(data_dir):
    importer = get_importer(data_dir / "IMMS" / "test_watersimms_txt.txt")
    data = importer.get_imms_avg_scan()
    assert data.ndim == 2 and data.shape[1] == 3
    np.testing.assert_array_equal(importer.get_imms_scan(1), data)
    assert len(importer.get_all_imms_scans()) == 1


def test_binary_pairs_and_case_insensitive_extension(tmp_path):
    path = tmp_path / "scan.BIN"
    np.array([[100.0, 2.0], [101.0, 3.0]]).tofile(path)
    importer = get_importer(path)
    np.testing.assert_array_equal(importer.get_avg_scan(), [[100, 2], [101, 3]])


def test_binary_odd_value_count_is_rejected(tmp_path):
    path = tmp_path / "scan.bin"
    np.array([1.0, 2.0, 3.0]).tofile(path)
    with pytest.raises(ValueError, match="pairs"):
        SingleScanImporter(path)


def test_two_column_file_is_not_ion_mobility(tmp_path):
    path = tmp_path / "scan.txt"
    np.savetxt(path, [[100, 2], [101, 3]])
    importer = SingleScanImporter(path)
    with pytest.raises(ValueError, match="ion-mobility"):
        importer.get_imms_avg_scan()


def test_cdms_defaults_missing_metadata_columns(tmp_path):
    path = tmp_path / "scan.csv"
    np.savetxt(path, [[100, 2], [101, 3]], delimiter=",")
    cdms = SingleScanImporter(path).get_cdms_data()
    np.testing.assert_array_equal(cdms[:, 2], [1, 1])
    np.testing.assert_array_equal(cdms[:, 3], [1, 1])
    np.testing.assert_array_equal(cdms[:, 4], [-1, -1])


def test_direct_single_scan_importer_rejects_unknown_extension(tmp_path):
    path = tmp_path / "scan.unknown"
    path.write_text("100 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported"):
        SingleScanImporter(path)
