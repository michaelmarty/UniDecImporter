# MassSpecImporter

MassSpecImporter reads open and vendor mass-spectrometry formats through one NumPy-based
API. It is independent of UniDec and runs on Windows, macOS, and Linux for open formats.

```python
from MassSpecImporter import get_importer

with get_importer("sample.mzML") as reader:
    data = reader.get_avg_scan(time_range=(1.0, 2.0))
```

Continue with [Getting started](getting-started.md), review the
[format matrix](formats.md), or open the [Python API](api.md).

