# Getting started

Install from PyPI:

```shell
python -m pip install mass-spec-importer
```

The factory accepts strings and path-like objects:

```python
from pathlib import Path
from MassSpecImporter import get_importer

with get_importer(Path("sample.mzML")) as reader:
    first_scan = reader.get_single_scan(reader.scans[0])
    average = reader.get_avg_scan(scan_range=(reader.scans[0], reader.scans[9]))
    print(reader.get_polarity(), average.shape)
```

Every ordinary spectrum is an `N x 2` NumPy array. CD-MS events are `N x 5`, and
ion-mobility data are `N x 3`. Scan ranges are inclusive.

Single-scan text files may begin with arbitrary header rows. TXT and DAT use whitespace;
CSV uses commas. NPZ files must contain an array named `data`. BIN files are interleaved
float64 m/z and intensity pairs.

Use the capability attributes before specialized operations:

```python
if reader.chrom_support:
    tic = reader.get_tic()
if reader.cdms_support:
    events = reader.get_cdms_data()
if reader.imms_support:
    mobility = reader.get_imms_avg_scan(mzbins=1)
```

