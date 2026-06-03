import csv
import importlib.util
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable
TEST_DATA_DIR = PROJECT_ROOT / "tests" / "data"


def require_runtime_dependencies():
    missing = [
        module
        for module in ("pytest", "pandas", "pyteomics", "Bio", "scipy", "lxml")
        if importlib.util.find_spec(module) is None
    ]
    if missing:
        pytest.skip("missing runtime dependencies: " + ", ".join(missing))


def require_external_analysis_fixtures():
    if not TEST_DATA_DIR.is_dir():
        pytest.skip(f"missing test data directory: {TEST_DATA_DIR}")


def run_cli(args):
    result = subprocess.run(
        [PYTHON, *args],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            "command failed\n"
            f"args: {args}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def extract_zip(zip_path, destination):
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(destination)


def read_tsv_rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            cleaned = {}
            for key, value in row.items():
                normalized_key = "" if key is None else key.strip()
                cleaned[normalized_key] = "" if value is None else value.strip()
            rows.append(cleaned)
        return rows


def build_single_spectrum_dir(tmp_path, spectrum_name="Rattus-TOF.csv"):
    extract_zip(PROJECT_ROOT / "Wiki_examples" / "Classify" / "Spectra.zip", tmp_path)
    source = tmp_path / "Spectra" / spectrum_name
    target_dir = tmp_path / "spectra"
    target_dir.mkdir()
    shutil.copy(source, target_dir / spectrum_name)
    return target_dir
