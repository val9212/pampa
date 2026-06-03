import pytest

from tests.cli_helpers import PROJECT_ROOT, TEST_DATA_DIR, extract_zip, read_tsv_rows, require_external_analysis_fixtures, require_runtime_dependencies, run_cli


@pytest.mark.integration
def test_craft_allpeptides_cli_with_spectra_filters_markers(tmp_path):
    require_runtime_dependencies()
    require_external_analysis_fixtures()

    spectra_dir = tmp_path / "allpeptides_spectra"
    extract_zip(TEST_DATA_DIR / "allpeptides_sp.zip", tmp_path)
    (tmp_path / "SPECTRA_SHEEP").rename(spectra_dir)

    output = tmp_path / "allpeptides_results.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--allpeptides",
            "-f",
            str(TEST_DATA_DIR / "allpeptides_s.fasta"),
            "-l",
            str(TEST_DATA_DIR / "limit.txt"),
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-c",
            str(PROJECT_ROOT / "config_mammals.json"),
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) > 20
    assert {row["Taxon name"] for row in rows} == {"Ovis aries"}
    assert any(row["Marker"] == "COL1A2-484-498" and row["Sequence"] == "GIPGEFGLPGPAGAR" for row in rows)
    assert any(row["Marker"] == "COL1A2-793-816" and row["Sequence"] == "GLPGVAGSVGEPGPLGIAGPPGAR" for row in rows)
    assert all(row["Digestion"] for row in rows)


@pytest.mark.integration
def test_craft_allpeptides_cli_without_spectra_generates_in_silico_table(tmp_path):
    require_runtime_dependencies()
    require_external_analysis_fixtures()

    output = tmp_path / "allpeptides_no_spectra.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--allpeptides",
            "-f",
            str(TEST_DATA_DIR / "allpeptides_s.fasta"),
            "-l",
            str(TEST_DATA_DIR / "limit.txt"),
            "-c",
            str(PROJECT_ROOT / "config_mammals.json"),
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) > 50
    assert {row["Taxon name"] for row in rows} == {"Ovis aries"}
    assert all(row["Comment"] == "In silico digestion." for row in rows)
