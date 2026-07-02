import shutil

import pytest

from tests.cli_helpers import PROJECT_ROOT, TEST_DATA_DIR, extract_zip, read_tsv_rows, require_external_analysis_fixtures, require_runtime_dependencies, run_cli


@pytest.mark.integration
def test_craft_fillin_cli_with_sequence_directory_completes_missing_fields(tmp_path):
    require_runtime_dependencies()
    require_external_analysis_fixtures()

    sequences_dir = tmp_path / "fillin_sequences"
    extract_zip(TEST_DATA_DIR / "fillin_s.zip", tmp_path)
    sequences_dir.mkdir()
    for fasta_file in ("fillin_s.fasta", "fillin_s2.fasta"):
        shutil.move(tmp_path / fasta_file, sequences_dir / fasta_file)

    output = tmp_path / "fillin_results.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--fillin",
            "-p",
            str(TEST_DATA_DIR / "fillin_pt.tsv"),
            "-d",
            str(sequences_dir),
            "-e",
            "0.1",
            "-t",
            str(TEST_DATA_DIR / "taxonomy.tsv"),
            "-c",
            str(PROJECT_ROOT / "config_mammals.json"),
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) == 5
    assert {row["Taxon name"] for row in rows} == {"Mus musculus"}
    supplemented_rows = [row for row in rows if row["Comment"].startswith("Supplement from")]

    assert len(supplemented_rows) == 5
    assert any(row["Sequence"] == "SGQPGPVGPAGVR" and row["PTM"] == "0H" for row in supplemented_rows)
    assert any(row["Sequence"] == "GLPGEFGLPGPAGPR" and row["Start"] == "580" for row in supplemented_rows)


@pytest.mark.integration
def test_craft_fillin_cli_with_fasta_file_completes_local_wiki_example(tmp_path):
    require_runtime_dependencies()

    output = tmp_path / "fillin_wiki.tsv"
    run_cli(
        [
            "pampa_craft.py",
            "--fillin",
            "-p",
            "Wiki_examples/table_mouse_ABC.tsv",
            "-f",
            "Wiki_examples/murine.fasta",
            "-e",
            "0.1",
            "-c",
            "config_mammals.json",
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    supplemented_rows = [row for row in rows if row["Comment"].startswith("Supplement from")]
    assert len(supplemented_rows) == 5
    assert any(row["Sequence"] == "SGQPGPVGPAGVR" and row["Start"] == "1074" for row in supplemented_rows)
    assert any(row["Sequence"] == "GATGLPGVAGAPGLPGPR" and row["Hel"] == "220" for row in supplemented_rows)
