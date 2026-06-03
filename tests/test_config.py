import json

from src import config


def test_config_accessors_read_expected_sections(tmp_path):
    config_file = tmp_path / "config.json"
    config_file.write_text(
        json.dumps(
            {
                "enzyme": "trypsin",
                "number_of_missed_cleavages": 2,
                "min_peptide_length": 8,
                "max_peptide_length": 40,
                "taxonomy_ranks": ["Kingdom", "Species"],
                "peptide_table_order": ["Marker", "Mass"],
                "marker_order": ["A", "B"],
                "min_number_of_peaks": 3,
                "min_proportion_of_peaks": 0.25,
                "taxonomy": "taxonomy.tsv",
                "peptide_table": ["table.tsv"],
                "substitution_matrices": ["x.csv"],
                "gamma_matrices": ["g.csv"],
                "conserved": "c.csv",
            }
        )
    )

    assert config.config_digestion(str(config_file)) == {
        "enzyme": "trypsin",
        "number_of_missed_cleavages": 2,
        "min_peptide_length": 8,
        "max_peptide_length": 40,
    }
    assert config.config_headers(str(config_file)) == [
        "Kingdom",
        "Species",
        "Marker",
        "Mass",
    ]
    assert config.config_markers(str(config_file)) == ["A", "B"]
    assert config.config_minimum_number_of_peaks(str(config_file)) == 3
    assert config.config_selection_peaks(str(config_file)) == 0.25
    assert config.config_taxonomy(str(config_file)) == "taxonomy.tsv"
    assert config.config_peptide_table(str(config_file)) == ["table.tsv"]
    assert config.config_matrices_and_co(str(config_file)) == (
        ["x.csv"],
        ["g.csv"],
        "c.csv",
    )
