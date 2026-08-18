# MassSpecImporter

MassSpecImporter is a standalone Python package for reading mass-spectrometry data into
NumPy arrays. It was extracted from UniDec, but has no runtime dependency on UniDec or
IsoDec.

The same high-level API covers single spectra, LC-MS, CD-MS, and ion-mobility MS:

```python
from MassSpecImporter import get_importer

with get_importer("run.mzML") as reader:
    spectrum = reader.get_avg_scan(time_range=(2.0, 2.5))
    tic = reader.get_tic()
```

## Installation

```shell
python -m pip install mass-spec-importer
```

For Thermo RAW support on Windows x86-64, install the optional bridge too:

```shell
python -m pip install "mass-spec-importer[thermo]"
```

Python 3.10–3.13 is supported on Windows, macOS, and Linux. Open formats work on all
three operating systems. The distributed package is pure Python; its open-format
dependencies provide normal platform wheels.

## Format support

| Format | Data types | Platforms | Reader dependency |
|---|---|---|---|
| mzML, indexed mzML.gz | LC-MS, MS/MS, CD-MS, IM-MS | Windows, macOS, Linux | pymzML |
| mzXML | LC-MS, MS/MS, CD-MS | Windows, macOS, Linux | pyteomics |
| TXT, DAT, CSV | Single scan, CD-MS, IM-MS | Windows, macOS, Linux | NumPy |
| NPZ, BIN | Single scan, CD-MS | Windows, macOS, Linux | NumPy |
| I2MS, DMT | CD-MS | Windows, macOS, Linux | Python sqlite3 |
| Thermo RAW | LC-MS, MS/MS, CD-MS | Windows x86-64 only | Thermo RawFileReader + pythonnet |
| Waters RAW directory | LC-MS, IM-MS | Windows x86-64 only | Waters MassLynx SDK |
| Agilent `.d` | LC-MS, MS/MS | Windows x86-64 only | Agilent MassHunter Data Access SDK + pythonnet |

Vendor SDK licenses prevent their binaries from being redistributed in the wheel or
source distribution. A local source checkout can use already-present licensed binaries.
For an installed package, set `THERMO_RAW_FILE_READER_DIR` to a licensed Thermo assembly
directory, `MASSLYNX_RAW_DLL` to a licensed Waters DLL, or `AGILENT_DA_SDK_DIR` to a
licensed Agilent Data Access assembly directory. Unsupported platforms and missing
runtimes raise `VendorReaderUnavailableError` with an actionable message.

## Core API

`get_importer(path, **options)` returns the appropriate reader. The common methods are:

- `get_single_scan(scan)` → `N x 2` array (`m/z`, intensity)
- `get_all_scans()` → list of `N x 2` arrays
- `get_avg_scan(scan_range=..., time_range=...)` → merged `N x 2` array
- `get_tic()` / `get_eic(mass, mz_tol, rt_range=None)` → chromatogram `N x 2` arrays
- `get_cdms_data()` → `N x 5` (`m/z`, intensity, scan, inverse injection time, time)
- `get_imms_scan(scan)` / `get_imms_avg_scan(...)` → `N x 3` (`m/z`, drift time, intensity)
- `get_polarity()`, `get_ms_order(scan)`, scan/time conversion helpers, and `close()`

Check `reader.chrom_support`, `reader.cdms_support`, and `reader.imms_support` before
calling specialized methods. Readers are context managers, so `with` is preferred.

For a two-column single spectrum:

```python
from MassSpecImporter import get_importer

reader = get_importer("spectrum.csv")
data = reader.get_avg_scan()
assert data.shape[1] == 2
```

For CD-MS:

```python
with get_importer("ions.dmt") as reader:
    events = reader.get_cdms_data()
    mz, intensity, scan, inverse_injection_time, time = events.T
```

## Development

Large test fixtures use Git LFS:

```shell
git clone https://github.com/michaelmarty/MassSpecImporter.git
cd MassSpecImporter
git lfs pull
python -m pip install -e ".[test]"
python -m pytest
```

The suite includes fast numerical/unit tests, cross-platform integration tests against
the bundled open formats, and separately marked Windows vendor tests:

```shell
python -m pytest -m "not integration and not vendor"
python -m pytest -m "integration and not vendor"
python -m pytest -m vendor
```

See the [documentation](https://michaelmarty.github.io/MassSpecImporter/),
[`PUBLISHING.md`](PUBLISHING.md), and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)
for full details.
