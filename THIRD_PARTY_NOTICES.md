# Third-party notices

MassSpecImporter itself is distributed under the BSD 3-Clause License.

## Thermo RawFileReader

Thermo support uses Thermo Fisher Scientific's RawFileReader .NET assemblies and
`pythonnet`. The RawFileReader license permits limited redistribution subject to
additional requirements and restrictions. To avoid imposing those terms on the
open-source distribution, the assemblies are **not included in wheels or source
distributions**. A copy of the vendor license is retained at
`MassSpecImporter/Thermo/RawFileReaderLicense.doc`.

Install a properly licensed RawFileReader locally and set
`THERMO_RAW_FILE_READER_DIR` to the directory containing its assemblies.

RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc.
All rights reserved.

## Waters MassLynx SDK

Waters' MassLynx SDK EULA prohibits redistribution. Consequently,
`MassLynxRaw.dll`, `cdt.dll`, and the import library are **not included in wheels
or source distributions**. The EULA supplied with the original SDK is retained at
`MassSpecImporter/Waters/Waters_MassLynxSDK_EULA.txt`.

Users must obtain their own licensed SDK and set `MASSLYNX_RAW_DLL` to the full
path of `MassLynxRaw.dll`.

## Agilent MassHunter Data Access SDK

Agilent support uses the proprietary MassHunter Data Access .NET assemblies. Those
assemblies are **not included in wheels or source distributions**. Users must obtain
and license them directly from Agilent, then set `AGILENT_DA_SDK_DIR` to their directory.

## Open-format dependencies

NumPy, h5py, pymzML, and pyteomics are installed as dependencies. Pythonnet is an
optional Thermo/Agilent/vendor dependency. Each remains governed by its respective license.
