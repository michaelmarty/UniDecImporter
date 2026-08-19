# Third-party notices

MassSpecImporter itself is distributed under the BSD 3-Clause License.
The BSD license does not apply to the bundled Thermo assemblies.

## Thermo RawFileReader

Thermo support uses Thermo Fisher Scientific's RawFileReader .NET assemblies and
`pythonnet`. The four RawFileReader assemblies used by this package are included in
wheels and source distributions under Thermo Fisher Scientific's separate proprietary
license. The authoritative license is included at
`MassSpecImporter/Thermo/RawFileReaderLicense.doc`; a readable notice and the terms that
apply to package users are in `THERMO_RAWFILEREADER_TERMS.md`.

Installing or using the bundled assemblies constitutes acceptance of those terms.
Among other restrictions, end users may not redistribute the Thermo assemblies,
commercially exploit them or products incorporating them without Thermo's prior written
consent, reverse engineer them, or use Thermo trademarks to imply endorsement. These
restrictions apply only to the Thermo assemblies, not to MassSpecImporter's BSD-licensed
Python source.

`THERMO_RAW_FILE_READER_DIR` may still be used to select another properly licensed
RawFileReader installation instead of the bundled copy.

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
