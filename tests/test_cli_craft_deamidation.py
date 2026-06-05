import pytest

from tests.cli_helpers import read_tsv_rows, require_runtime_dependencies, run_cli, TEST_DATA_DIR


@pytest.mark.integration
def test_craft_deamidation_cli_adds_modified_rows(tmp_path):
    require_runtime_dependencies()

    output = tmp_path / "deamidation_results.tsv"
    limit_file = tmp_path / "limit.txt"
    limit_file.write_text("Deamidation = G\n")

    run_cli(
        [
            "pampa_craft.py",
            "--deamidation",
            "-p",
            str(TEST_DATA_DIR / "deamidation.tsv"),
            "-l",
            str(limit_file),
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) == 4
    assert any(row["PTM"] == "4H1D" for row in rows)
    assert any(row["PTM"] == "5H1D" for row in rows)
    assert any(row["Mass"] == "3000.484016" for row in rows)
    assert any(row["Mass"] == "2984.49597168482" for row in rows)
