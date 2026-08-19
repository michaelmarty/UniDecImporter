# Formats and platforms

## Cross-platform formats

- mzML and indexed mzML.gz: LC-MS, MS/MS, CD-MS, and ion mobility through pymzML.
- mzXML: LC-MS, MS/MS, and CD-MS through pyteomics.
- TXT, DAT, and CSV: single spectra plus three-column ion-mobility and multi-column
  CD-MS data.
- NPZ: the `data` member may contain two or more columns.
- BIN: interleaved float64 m/z and intensity pairs.
- I2MS and DMT: SQLite ion-event databases.

These readers are tested on Windows, macOS, and Linux. Open-format modules are imported
without loading any vendor runtime.

## Windows vendor formats

Thermo RAW and Waters RAW require 64-bit x86 Windows. The factory distinguishes a Thermo
RAW file from a Waters RAW directory.

The Windows-only Thermo RawFileReader assemblies are bundled as data files under their
separate proprietary license. Waters and Agilent binaries are not redistributed; obtain
and license those SDKs from their vendors. The following variables override the bundled
Thermo location or configure a locally installed vendor SDK:

```powershell
$env:THERMO_RAW_FILE_READER_DIR = "C:\path\to\RawFileReader"
$env:MASSLYNX_RAW_DLL = "C:\path\to\MassLynxRaw.dll"
$env:AGILENT_DA_SDK_DIR = "C:\path\to\MassHunter\DataAnalysis"
```

Thermo and Agilent also require `pythonnet`; install it with
`python -m pip install "mass-spec-importer[vendor]"` on Windows AMD64.
Waters companion DLLs should be placed beside `MassLynxRaw.dll` or otherwise be visible
to the Windows loader.

Agilent `.d` reads spectra and chromatograms through the MassHunter Data Access SDK.
