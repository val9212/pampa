import importlib
import json
import sys
import types

import pytest

from src import message


def import_params_checker():
    sys.modules.pop("src.params_checker", None)
    sys.modules["src.report"] = types.SimpleNamespace(
        create_report_header=lambda *args, **kwargs: None
    )
    return importlib.import_module("src.params_checker")


def test_logger_and_outputdir_configuration_creates_directory(tmp_path):
    params_checker = import_params_checker()

    output = tmp_path / "nested" / "results"
    output_dir, output_file, report_file, detail_file, output_json = (
        params_checker.logger_and_outputdir_configuration(str(output), "cmd")
    )

    assert (tmp_path / "nested").is_dir()
    assert output_dir == str(tmp_path / "nested")
    assert output_file == "results.tsv"
    assert report_file == "report_results.txt"
    assert detail_file == "detail_results.tsv"
    assert output_json == "results.json"


def test_check_and_update_parameters_classify_resets_invalid_neighbour(tmp_path):
    params_checker = import_params_checker()

    spectra_dir = tmp_path / "spectra"
    spectra_dir.mkdir()
    peptide_table = tmp_path / "table.tsv"
    peptide_table.write_text("dummy")
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({"taxonomy": "taxonomy.tsv", "peptide_table": []}))

    message.configure(str(tmp_path))
    result = params_checker.check_and_update_parameters_classify(
        spectra=str(spectra_dir),
        taxonomy=None,
        peptide_table=[str(peptide_table)],
        fasta=None,
        fasta_dir=None,
        limit=None,
        deamidation=False,
        error=0.1,
        neighbour=150,
        allpeptides=False,
        mammals=False,
        placentals=False,
        birds=False,
        config=str(config_file),
    )

    assert result[8] == 100
    assert "Parameter -n (neighbouring): value is 100" in (
        tmp_path / "warning.log"
    ).read_text()


def test_check_and_update_parameters_classify_rejects_missing_marker_source(tmp_path):
    params_checker = import_params_checker()

    spectra_dir = tmp_path / "spectra"
    spectra_dir.mkdir()
    config_file = tmp_path / "config.json"
    config_file.write_text(json.dumps({"taxonomy": "taxonomy.tsv", "peptide_table": []}))

    message.configure(str(tmp_path))
    with pytest.raises(message.InputError):
        params_checker.check_and_update_parameters_classify(
            spectra=str(spectra_dir),
            taxonomy=None,
            peptide_table=None,
            fasta=None,
            fasta_dir=None,
            limit=None,
            deamidation=False,
            error=0.1,
            neighbour=100,
            allpeptides=False,
            mammals=False,
            placentals=False,
            birds=False,
            config=str(config_file),
        )
