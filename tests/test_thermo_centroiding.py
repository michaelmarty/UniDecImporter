import platform
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest


pytestmark = pytest.mark.skipif(
    platform.system() != "Windows",
    reason="Thermo RawFileReader requires Windows",
)


def test_average_spectrum_can_force_thermo_centroid_stream(monkeypatch):
    raw_module = pytest.importorskip("UniDecImporter.Thermo.RawFileReader")
    options = SimpleNamespace(ToleranceUnits=None, Tolerance=None)
    centroid = SimpleNamespace(
        Masses=[100.0, 101.0],
        Intensities=[10.0, 20.0],
        Charges=[1, 1],
        Length=2,
    )
    average = SimpleNamespace(CentroidScan=centroid)
    calls = []

    class FakeExtensions:
        @staticmethod
        def DefaultMassOptions(source):
            return options

        @staticmethod
        def BuildFilterHelper(source, name):
            return SimpleNamespace(Filter=name)

        @staticmethod
        def AverageScansInScanRange(source, first, last, scan_filter, mass_options):
            calls.append((first, last, scan_filter, mass_options))
            return average

    monkeypatch.setattr(raw_module, "Extensions", FakeExtensions)
    monkeypatch.setattr(raw_module, "ToleranceUnits", SimpleNamespace(ppm="ppm"))
    reader = raw_module.RawFileReader.__new__(raw_module.RawFileReader)
    reader.source = SimpleNamespace(
        GetScanStatsForScanNumber=lambda scan: SimpleNamespace(IsCentroidScan=False),
    )

    result = reader.GetAverageSpectrum([2, 4], centroid=True)

    np.testing.assert_array_equal(result, [[100.0, 10.0], [101.0, 20.0]])
    assert calls == [(2, 4, "FTMS", options)]
    assert options.ToleranceUnits == "ppm"
    assert options.Tolerance == 5.0


def test_grab_all_centroid_data_uses_averaged_centroid_stream():
    thermo_module = pytest.importorskip("UniDecImporter.Thermo.Thermo")
    expected = np.array([[5.0, 1.0], [100.0, 10.0], [101.0, 20.0]])
    importer = thermo_module.ThermoImporter.__new__(thermo_module.ThermoImporter)
    importer.scan_range = np.array([2, 4])
    importer.msrun = Mock()
    importer.msrun.GetAverageSpectrum.return_value = expected

    result = importer.grab_all_centroid_dat()

    np.testing.assert_array_equal(result, expected[1:])
    args, kwargs = importer.msrun.GetAverageSpectrum.call_args
    np.testing.assert_array_equal(args[0], importer.scan_range)
    assert kwargs == {"centroid": True}
    importer.msrun.GetCentroidSpectrum.assert_not_called()
