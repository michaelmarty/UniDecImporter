# Python API

## Factory and errors

`MassSpecImporter.get_importer(file_path, **kwargs)` delegates to
`ImporterFactory.create_importer`. Unknown extensions raise `UnsupportedFormatError`.
Unavailable vendor platforms or runtimes raise `VendorReaderUnavailableError`.

`MassSpecImporter.get_polarity(path)` is a convenience function that returns a polarity
string and preserves the historical positive fallback on read failure.

## Common reader state

- `scans`: ordered scan identifiers.
- `times`: retention times in minutes.
- `levels`: MS order where the format provides it.
- `scan_range`: inclusive minimum and maximum scan identifiers.
- `chrom_support`, `cdms_support`, `imms_support`: capability flags.

## Common methods

- `get_single_scan(scan)` returns a two-column spectrum.
- `get_all_scans()` returns spectra in scan order.
- `get_avg_scan(scan_range=None, time_range=None)` merges an inclusive range.
- `get_scan_time(scan)`, `get_time_scan(time)`, `get_scans_from_times(range)`, and
  `get_times_from_scans(range)` translate scan and retention-time coordinates.
- `get_tic()` and `get_eic(mass, mz_tol, rt_range=None)` return time/intensity arrays.
- `get_cdms_data()` returns m/z, intensity, scan, inverse injection time, and time.
- `get_imms_scan(scan)` and `get_imms_avg_scan(...)` return m/z, drift time, intensity.
- `check_centroided()` applies the historical lag-one autocorrelation heuristic.
- `close()` releases native or file resources. All readers support context-manager use.

## FileParser

`MassSpecImporter.FileParser` contains batch slicing helpers for writing averaged time or
scan windows to text or UniDec-compatible HDF5 datasets: `parse`, `parse_multiple`,
`extract`, `extract_scans`, `extract_timepoints`, and `extract_scans_multiple_files`.

