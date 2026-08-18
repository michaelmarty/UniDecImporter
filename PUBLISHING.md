# Building and publishing MassSpecImporter

MassSpecImporter uses setuptools and produces a universal pure-Python wheel plus a source
distribution. Vendor DLLs and the large Git LFS test corpus are excluded from both.

## Prepare a release

1. Update `MassSpecImporter/_version.py` and `CITATION.cff` to the same version.
2. Update the README, documentation, release notes, and third-party notices.
3. Run `git lfs pull` and the complete test suite on licensed development machines.
4. Build locally and inspect the artifacts for accidental DLL or test-data inclusion.
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

Use `wheel-test/bin/python` on macOS and Linux. Verify that neither archive contains
`TestData`, `*.dll`, or `*.lib` before publishing.

The workflow creates a GitHub release on manual dispatch and can publish to PyPI through
trusted publishing. Configure a PyPI environment named `pypi`; no long-lived token is
required.

