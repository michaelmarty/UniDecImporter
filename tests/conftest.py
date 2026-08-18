from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "TestData"


@pytest.fixture(scope="session")
def data_dir():
    if not DATA.exists():
        pytest.skip("bundled integration data are not available")
    return DATA


def require_real_data(path: Path, minimum_size=1024):
    """Skip cleanly when a Git LFS object was not downloaded."""
    if not path.exists() or (path.is_file() and path.stat().st_size < minimum_size):
        pytest.skip(f"test data unavailable (run git lfs pull): {path}")
    return path

