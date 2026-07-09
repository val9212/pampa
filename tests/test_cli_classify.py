import shutil

import pytest

from tests.cli_helpers import PROJECT_ROOT, build_single_spectrum_dir, read_tsv_rows, require_runtime_dependencies, run_cli


@pytest.mark.integration
def test_classify_cli_matches_wiki_example(tmp_path):
    require_runtime_dependencies()

    spectra_dir = tmp_path / "spectra"
    from tests.cli_helpers import extract_zip

    extract_zip(PROJECT_ROOT / "Wiki_examples" / "Classify" / "Spectra.zip", tmp_path)
    (tmp_path / "Spectra").rename(spectra_dir)

    output = tmp_path / "classify_results.tsv"
    run_cli(
        [
            "pampa_classify.py",
            "--mammals",
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-o",
            str(output),
        ]
    )

    assert output.is_file()
    assert (tmp_path / "detail_classify_results.tsv").is_file()
    assert (tmp_path / "report_classify_results.txt").is_file()

    rows = read_tsv_rows(output)
    assert len(rows) == 6

    by_spectrum = {row["Spectrum"]: row for row in rows}
    assert by_spectrum["Castor-TOF.csv"]["Assignment"] == "51338 [Castor canadensis]"
    assert by_spectrum["Castor-TOF.csv"]["Maximal clade"] == "29132 [Castoridae]"
    assert by_spectrum["Horse-TOF.csv"]["Assignment"] == "9789 [Equus]"
    assert "Equus caballus" in by_spectrum["Horse-TOF.csv"]["Species"]
    assert by_spectrum["Whale-TOF.csv"]["Assignment"] == "9766 [Balaenoptera]"


@pytest.mark.integration
def test_classify_cli_accepts_custom_peptide_table_and_taxonomy(tmp_path):
    require_runtime_dependencies()

    spectra_dir = build_single_spectrum_dir(tmp_path)
    output = tmp_path / "classify_pt.tsv"
    run_cli(
        [
            "pampa_classify.py",
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-p",
            "Peptide_tables/table_mammals.tsv",
            "-t",
            "Taxonomy/taxonomy_mammals.tsv",
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) == 1
    assert rows[0]["Spectrum"] == "Rattus-TOF.csv"
    assert rows[0]["Assignment"] == "10116 [Rattus norvegicus]"
    assert rows[0]["Rank"] == "species"


@pytest.mark.integration
def test_classify_cli_accepts_fasta_input(tmp_path):
    require_runtime_dependencies()

    spectra_dir = build_single_spectrum_dir(tmp_path)
    output = tmp_path / "classify_fasta.tsv"
    run_cli(
        [
            "pampa_classify.py",
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-f",
            "Wiki_examples/murine.fasta",
            "-t",
            "Taxonomy/taxonomy_mammals.tsv",
            "-c",
            "config_mammals.json",
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) == 1
    assert rows[0]["Spectrum"] == "Rattus-TOF.csv"
    assert rows[0]["Assignment"] == "10090 [Mus musculus]"
    assert rows[0]["#peaks"] == "35"


@pytest.mark.integration
def test_classify_cli_accepts_fasta_directory_input(tmp_path):
    require_runtime_dependencies()

    spectra_dir = build_single_spectrum_dir(tmp_path)
    fasta_dir = tmp_path / "murine_dir"
    fasta_dir.mkdir()
    shutil.copy(PROJECT_ROOT / "Wiki_examples" / "murine.fasta", fasta_dir / "murine.fasta")

    output = tmp_path / "classify_dir.tsv"
    run_cli(
        [
            "pampa_classify.py",
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-d",
            str(fasta_dir),
            "-t",
            "Taxonomy/taxonomy_mammals.tsv",
            "-c",
            "config_mammals.json",
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) == 1
    assert rows[0]["Spectrum"] == "Rattus-TOF.csv"
    assert rows[0]["Assignment"] == "10090 [Mus musculus]"


@pytest.mark.integration
def test_classify_cli_accepts_suboptimal_and_isotope_flags(tmp_path):
    require_runtime_dependencies()

    spectra_dir = build_single_spectrum_dir(tmp_path)
    output = tmp_path / "classify_pt_options.tsv"
    run_cli(
        [
            "pampa_classify.py",
            "-s",
            str(spectra_dir),
            "-e",
            "0.1",
            "-p",
            "Peptide_tables/table_mammals.tsv",
            "-t",
            "Taxonomy/taxonomy_mammals.tsv",
            "-n",
            "80",
            "-a",
            "-i",
            "-o",
            str(output),
        ]
    )

    rows = read_tsv_rows(output)
    assert len(rows) >= 1
    assert rows[0]["Spectrum"] == "Rattus-TOF.csv"
    assert (tmp_path / "detail_classify_pt_options.tsv").is_file()
