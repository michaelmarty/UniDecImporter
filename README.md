# UniDecImporter

UniDecImporter is a standalone Python package for reading mass-spectrometry data into
NumPy arrays. It was extracted from UniDec, but has no runtime dependency on UniDec or
IsoDec.

The same high-level API covers single spectra, LC-MS, CD-MS, and ion-mobility MS:

```python
from UniDecImporter import get_importer

with get_importer("run.mzML") as reader:
    spectrum = reader.get_avg_scan(time_range=(2.0, 2.5))
    tic = reader.get_tic()
```

## Installation

```shell
python -m pip install UniDecImporter
```

For Thermo RAW support on Windows x86-64, install the optional bridge too:

```shell
python -m pip install "UniDecImporter[thermo]"
```

Python 3.10–3.13 is supported on Windows, macOS, and Linux. Open formats work on all
three operating systems. The package includes Windows-only Thermo .NET assemblies as
data files, but importing and using open-format readers does not load them.

> **Thermo proprietary software:** Installing or using the bundled Thermo RawFileReader
> assemblies means you accept Thermo's separate license, included in the distribution.
> End users may not redistribute those assemblies. Commercial exploitation requires
> Thermo's prior written consent. The BSD license covers this project's Python code, not
> the Thermo binaries. See [`THERMO_RAWFILEREADER_TERMS.md`](THERMO_RAWFILEREADER_TERMS.md).
>
> RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All
> rights reserved.

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

The four Thermo RawFileReader assemblies are bundled under their separate proprietary
terms. Waters and Agilent SDK licenses do not permit this project to redistribute those
binaries, so users must obtain them separately. Set `THERMO_RAW_FILE_READER_DIR` to
override the bundled Thermo assembly directory, `MASSLYNX_RAW_DLL` to a licensed Waters
DLL, or `AGILENT_DA_SDK_DIR` to a licensed Agilent Data Access assembly directory.
Unsupported platforms and missing runtimes raise `VendorReaderUnavailableError` with an
actionable message.

## Core API

`get_importer(file_path, **options)` inspects the path and returns the appropriate
reader. `file_path` may be a string or path-like object; `options` are passed to the
format-specific reader. Unknown extensions raise `UnsupportedFormatError`, while an
unavailable proprietary SDK or platform raises `VendorReaderUnavailableError`.

All readers inherit from `UniDecImporter.Importer.Importer`. They expose `scans`
(ordered scan identifiers), `times` (retention times in minutes), `levels` (MS order
when available), and an inclusive `scan_range`. Ordinary spectra and chromatograms are
NumPy `N x 2` arrays. Specialized return shapes are noted below.

### Spectrum access

| Method | Description |
|---|---|
| `get_single_scan(scan)` | Return one spectrum as `N x 2` (`m/z`, intensity). |
| `get_all_scans()` | Load and return all spectra as a list in scan order. |
| `get_avg_scan(scan_range=None, time_range=None)` | Merge spectra over an inclusive scan range or a retention-time range. Reader implementations may use a streaming path for large files. |
| `avg_fast(scan_range=None, time_range=None)` | Merge spectra through the shared cached-data implementation used by format readers. Most applications should call `get_avg_scan`. |
| `get_mz_localmax(mz, mz_tol)` | Find each scan's local peak near `mz`; `mz_tol` is in ppm. Return (`m/z`, intensity) rows and omit scans with no peak. |

Supplying `time_range` selects the nearest endpoint scans and takes precedence over a
simultaneously supplied `scan_range`. With neither argument, averaging uses the full
reader range.

### Scan metadata and coordinate conversion

| Method | Description |
|---|---|
| `get_polarity(scan=None)` | Return `"Positive"`, `"Negative"`, or a reader-specific unknown value. |
| `get_ms_order(scan=1)` | Return the scan's MS level; readers without level metadata default to `1`. |
| `get_max_scan()` | Return the final scan identifier. |
| `get_max_time()` | Return the final retention time in minutes. |
| `get_scan_index(scan)` | Convert a scan identifier to its array index. Identifiers outside the available range clamp to the first or last index. |
| `get_scan_time(scan)` | Return the retention time for a scan identifier. |
| `get_time_scan(time)` | Return the scan identifier nearest to a retention time in minutes. |
| `get_scans_from_times(time_range)` | Convert two retention times to an inclusive pair of scan identifiers. |
| `get_times_from_scans(scan_range)` | Return `[start, midpoint, end]` retention times for an inclusive scan range. |
| `scan_range_from_inputs(scan_range=None, time_range=None)` | Resolve either selection form and clamp it to the available scans. This is primarily useful to reader implementers. |
| `check_centroided()` | Estimate centroid/profile status from the first selected scan using a lag-one autocorrelation heuristic and update `reader.centroided`. |

### Chromatograms

These methods require `reader.chrom_support`:

| Method | Description |
|---|---|
| `get_tic()` | Return the total-ion chromatogram as `N x 2` (retention time, intensity). |
| `index_scans(min_mz, bin_width)` | Build or replace the in-memory MS1 peak index used for extracted-ion chromatograms. |
| `get_eic(mass, mz_tol, rt_range=None)` | Return an extracted-ion chromatogram around `mass`; `mz_tol` is an absolute m/z tolerance and `rt_range` is optional. The default index is built automatically on first use. |

### CD-MS and ion mobility

Check `reader.cdms_support` and `reader.imms_support` before calling these methods:

| Method | Description |
|---|---|
| `get_cdms_data(scan_range=None)` | Return an `N x 5` array: m/z, intensity, scan, inverse injection time, and retention time. Some readers override the base metadata defaults. |
| `get_imms_scan(scan)` | Return one ion-mobility scan as `N x 3` (m/z, drift time, intensity). |
| `get_all_imms_scans()` | Load all ion-mobility scans in scan order. |
| `get_imms_avg_scan(scan_range=None, time_range=None, mzbins=1)` | Merge selected mobility scans into an `N x 3` array. `mzbins` sets the linear m/z bin width; a false-like value requests automatic spacing. |

Calling an unsupported specialized method raises `NotImplementedError`.

### Construction and resource management

`Importer.__init__(file_path, **kwargs)` initializes the shared state, but applications
should construct readers through `get_importer`. `close()` releases open files, database
connections, or vendor handles; the base implementation is a no-op. `__enter__()` and
`__exit__()` provide context-manager support, call `close()` on exit, and do not suppress
exceptions. Prefer `with` whenever practical:

```python
with get_importer("run.mzML") as reader:
    data = reader.get_avg_scan()
```

For a two-column single spectrum:

```python
from UniDecImporter import get_importer

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
git clone https://github.com/michaelmarty/UniDecImporter.git
cd UniDecImporter
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

See the [documentation](https://michaelmarty.github.io/UniDecImporter/),
[`PUBLISHING.md`](PUBLISHING.md), and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)
for full details.
