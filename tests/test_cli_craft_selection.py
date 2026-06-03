import pytest

from tests.cli_helpers import build_single_spectrum_dir, read_tsv_rows, require_runtime_dependencies, run_cli


@pytest.mark.integration
def test_craft_selection_cli_filters_marker_table_against_spectra(tmp_path):
    require_runtime_dependencies()

    spectra_dir = build_single_spectrum_dir(tmp_path)
    output = tmp_path / "selection.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--selection",
            "-p",
            "Peptide_tables/table_mammals.tsv",
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-c",
            "config_mammals.json",
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) > 0
    assert all(row["Status"] == "MS" for row in rows)
    assert any("Rattus-TOF.csv" in row["Comment"] for row in rows)
