# Development

Clone the repository, download the Git LFS fixtures, and install an editable environment:

```shell
git clone https://github.com/michaelmarty/MassSpecImporter.git
cd MassSpecImporter
git lfs pull
python -m pip install -e ".[test]"
python -m pytest
```

On a licensed Windows vendor-reader workstation, install
`python -m pip install -e ".[test,vendor]"` before running the vendor-marked tests.

The tests are grouped as follows:

- unmarked tests cover numerical helpers, factory dispatch, batch parsing, error paths,
  and generated small files;
- `integration` tests read the full bundled data from every open format;
- `vendor` tests require Windows x86-64 and a locally licensed SDK.

Build strict documentation with:

```shell
python -m pip install -e ".[docs]"
python -m mkdocs build --strict
```

Build distributions with:

```shell
python -m pip install -e ".[release]"
python -m build
python -m twine check dist/*
```

The distribution is deliberately pure Python and excludes test data and all vendor DLLs.
See `PUBLISHING.md` for the release checklist.
