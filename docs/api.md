# Python API

## Factory and errors

`UniDecImporter.get_importer(file_path, **kwargs)` delegates to
`ImporterFactory.create_importer`. Unknown extensions raise `UnsupportedFormatError`.
Unavailable vendor platforms or runtimes raise `VendorReaderUnavailableError`.

`UniDecImporter.get_polarity(path)` is a convenience function that returns a polarity
string and preserves the historical positive fallback on read failure.

## Common reader state

- `scans`: ordered scan identifiers.
- `times`: retention times in minutes.
- `levels`: MS order where the format provides it.
- `scan_range`: inclusive minimum and maximum scan identifiers.
- `chrom_support`, `cdms_support`, `imms_support`: capability flags.

## Base `Importer` methods

### Spectra

- `get_single_scan(scan)` returns one `N x 2` m/z-intensity array.
- `get_all_scans()` returns a list of spectra in scan order.
- `get_avg_scan(scan_range=None, time_range=None)` merges an inclusive range. A supplied
  time range takes precedence and is converted to its nearest endpoint scans.
- `avg_fast(scan_range=None, time_range=None)` is the shared cached-data averaging
  implementation used by format readers.
- `get_mz_localmax(mz, mz_tol)` returns local (`m/z`, intensity) peaks across scans,
  using a tolerance in ppm.

### Metadata and range conversion

- `get_polarity(scan=None)` returns the ion polarity.
- `get_ms_order(scan=1)` returns the MS level, defaulting to MS1 when unavailable.
- `get_max_scan()` and `get_max_time()` return the final scan identifier and retention
  time in minutes.
- `get_scan_index(scan)` maps an identifier to an array index and clamps values outside
  the available range; `get_scan_time(scan)` returns the corresponding time.
- `get_time_scan(time)` returns the scan nearest to a time in minutes.
- `get_scans_from_times(time_range)` converts times to inclusive scan endpoints.
- `get_times_from_scans(scan_range)` returns start, midpoint, and end times.
- `scan_range_from_inputs(scan_range=None, time_range=None)` resolves and clamps either
  selection form for reader implementations.
- `check_centroided()` updates and returns `centroided` using the historical lag-one
  autocorrelation heuristic.

### Chromatograms

- `get_tic()` returns `N x 2` retention-time/intensity data.
- `index_scans(min_mz, bin_width)` builds the in-memory MS1 peak index used by EICs.
- `get_eic(mass, mz_tol, rt_range=None)` returns an extracted-ion chromatogram using an
  absolute m/z tolerance and builds a default index on first use.

These operations require `chrom_support`.

### CD-MS and ion mobility

- `get_cdms_data(scan_range=None)` returns `N x 5`: m/z, intensity, scan, inverse
  injection time, and retention time.
- `get_imms_scan(scan)` returns one `N x 3` m/z, drift-time, intensity array.
- `get_all_imms_scans()` returns every ion-mobility scan in order.
- `get_imms_avg_scan(scan_range=None, time_range=None, mzbins=1)` merges mobility scans;
  `mzbins` is the linear m/z bin width.

Check `cdms_support` or `imms_support` first. Unsupported operations raise
`NotImplementedError`.

### Lifetime

`Importer.__init__(file_path, **kwargs)` initializes common state, though callers
normally use the factory. `close()` releases resources. `__enter__()` returns the reader and
`__exit__()` closes it without suppressing exceptions, so all readers support `with`.

## FileParser

`UniDecImporter.FileParser` contains batch slicing helpers for writing averaged time or
scan windows to text or UniDec-compatible HDF5 datasets: `parse`, `parse_multiple`,
`extract`, `extract_scans`, `extract_timepoints`, and `extract_scans_multiple_files`.

