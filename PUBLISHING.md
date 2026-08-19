# Building and publishing MassSpecImporter

MassSpecImporter uses setuptools and produces a universal wheel plus a source
distribution. Both intentionally contain four Windows-only Thermo RawFileReader DLLs
under Thermo's separate proprietary license. Waters and Agilent binaries and the large
Git LFS test corpus are excluded.

Publishing the Thermo assemblies invokes the vendor license's distributor obligations,
including its indemnification requirement, required copyright notice, end-user
non-redistribution agreement, and restriction on commercial exploitation. Review the
authoritative `MassSpecImporter/Thermo/RawFileReaderLicense.doc` before every release;
`THERMO_RAWFILEREADER_TERMS.md` is only a readable summary and end-user notice.

## Prepare a release

1. Update `MassSpecImporter/_version.py` and `CITATION.cff` to the same version.
2. Update the README, documentation, release notes, and third-party notices.
3. Run `git lfs pull` and the complete test suite on licensed development machines.
4. Build locally and confirm that the four approved Thermo DLLs are present and that no
   other vendor binary or test data is included.
5. Push the release commit and run the **Build and publish** workflow.

## Local verification

```shell
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
```

Test the wheel from outside the checkout:

```shell
python -m venv wheel-test
wheel-test/Scripts/python -m pip install dist/mass_spec_importer-*.whl
cd wheel-test
Scripts/python -c "import MassSpecImporter; print(MassSpecImporter.__version__)"
```

Use `wheel-test/bin/python` on macOS and Linux. Verify that both archives contain exactly
these vendor binaries and no `TestData`, `*.lib`, or `*.tlb` files:

- `ThermoFisher.CommonCore.BackgroundSubtraction.dll`
- `ThermoFisher.CommonCore.Data.dll`
- `ThermoFisher.CommonCore.MassPrecisionEstimator.dll`
- `ThermoFisher.CommonCore.RawFileReader.dll`

Publishing must also retain `THERMO_RAWFILEREADER_TERMS.md` and
`MassSpecImporter/Thermo/RawFileReaderLicense.doc` in the artifacts.

The workflow creates a GitHub release on manual dispatch and can publish to PyPI through
trusted publishing. Configure a PyPI environment named `pypi`; no long-lived token is
required.
