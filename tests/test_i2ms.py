import numpy as np
import pytest

from UniDecImporter import get_importer
from UniDecImporter.I2MS.I2MS import I2MSImporter


@pytest.fixture
def i2ms(data_dir):
    importer = get_importer(data_dir / "CDMS" / "test_dmt_cdms.dmt")
    yield importer
    importer.close()


@pytest.mark.integration
def test_i2ms_metadata_and_scans(i2ms):
    assert isinstance(i2ms, I2MSImporter)
    assert i2ms.get_scan_range() == i2ms.scan_range
    assert i2ms.scan_range[1] > i2ms.scan_range[0]
    scans = i2ms.get_all_scans()
    assert len(scans) == len(i2ms.scans)
    assert all(scan.shape[1] == 2 for scan in scans)
    assert i2ms.get_single_scan(i2ms.scans[0]).shape[1] == 2
    assert sum(len(scan) for scan in i2ms.get_all_scans(threshold=0)) <= sum(map(len, scans))


@pytest.mark.integration
def test_i2ms_cdms_filter_and_average(i2ms):
    cdms = i2ms.get_cdms_data()
    assert cdms.shape[1] == 5 and len(cdms) > 100
    subset = i2ms.get_cdms_data_by_scans(i2ms.scan_range[:])
    np.testing.assert_array_equal(subset, cdms)
    average = i2ms.get_avg_scan(bins=2)
    assert average.shape[1] == 2
    assert average[:, 1].sum() == pytest.approx(cdms[:, 1].sum())
    with pytest.raises(ValueError):
        i2ms.get_avg_scan(bins=0)


@pytest.mark.integration
def test_i2ms_close_is_idempotent(i2ms):
    i2ms.close()
    i2ms.close()
