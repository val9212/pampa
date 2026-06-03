import shutil

import pytest

from tests.cli_helpers import PROJECT_ROOT, TEST_DATA_DIR, extract_zip, read_tsv_rows, require_external_analysis_fixtures, require_runtime_dependencies, run_cli


@pytest.mark.integration
def test_craft_homology_cli_with_fasta_generates_marker_table(tmp_path):
    require_runtime_dependencies()
    require_external_analysis_fixtures()

    output = tmp_path / "homology_results.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--homology",
            "-p",
            str(TEST_DATA_DIR / "homology_pt.tsv"),
            "-f",
            str(TEST_DATA_DIR / "homology_s.fasta"),
            "-t",
            str(TEST_DATA_DIR / "taxonomy.tsv"),
            "-c",
            str(PROJECT_ROOT / "config_mammals.json"),
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) >= 8
    assert {row["Taxon name"] for row in rows} == {"Lutra lutra"}
    assert {row["Marker"] for row in rows} >= {"COL1A1-508/P1", "COL1A2-978/A", "COL1A2-757/G"}
    assert any(row["Sequence"] == "GVQGPPGPAGPR" and row["PTM"] == "1H" for row in rows)


@pytest.mark.integration
def test_craft_homology_cli_with_sequence_directory_generates_marker_table(tmp_path):
    require_runtime_dependencies()
    require_external_analysis_fixtures()

    sequences_dir = tmp_path / "homology_sequences"
    extract_zip(TEST_DATA_DIR / "homology_s.zip", tmp_path)
    sequences_dir.mkdir()
    shutil.move(tmp_path / "homology_s.fasta", sequences_dir / "homology_s.fasta")

    output = tmp_path / "homology_dir_results.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--homology",
            "-p",
            str(TEST_DATA_DIR / "homology_pt.tsv"),
            "-d",
            str(sequences_dir),
            "-t",
            str(TEST_DATA_DIR / "taxonomy.tsv"),
            "-c",
            str(PROJECT_ROOT / "config_mammals.json"),
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) >= 8
    assert {row["Taxon name"] for row in rows} == {"Lutra lutra"}
    assert any(row["Sequence"] == "TGHPGTVGPAGIR" and row["PTM"] == "0H" for row in rows)
