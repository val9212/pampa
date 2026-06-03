import importlib
import sys
import types


def import_limit_module():
    sys.modules.pop("src.limit", None)
    sys.modules["src.taxonomy"] = types.SimpleNamespace(
        search_taxid_from_taxon_name=lambda name, taxonomy: name
    )
    return importlib.import_module("src.limit")


def test_parse_sequence_info_supports_filename_and_key_values():
    limit = import_limit_module()

    assert limit.parse_sequence_info("file_a.fasta", "limits.txt") == {
        "FileName": "file_a.fasta"
    }
    assert limit.parse_sequence_info(
        "OS = Homo sapiens, Pan troglodytes OX = 9606, 9598 PTM = 1H, 2H",
        "limits.txt",
    ) == {
        "OS": {"Homosapiens", "Pantroglodytes"},
        "OX": {"9606", "9598"},
        "PTM": {"1H", "2H"},
    }


def test_parse_limits_and_deamidation_selection(tmp_path):
    limit = import_limit_module()
    limit_file = tmp_path / "limits.txt"
    limit_file.write_text("Seq1.fasta\nDeamidation = A, B\nOS = Homo sapiens\n")

    parsed = limit.parse_limits(str(limit_file))

    assert parsed == [
        {"FileName": "Seq1.fasta"},
        {"Deamidation": {"A", "B"}},
        {"OS": {"Homosapiens"}},
    ]
    assert limit.deamidated_codes(parsed, False) == set()
    assert limit.deamidated_codes(parsed, True) == {"A", "B"}
    assert limit.deamidated_codes([{"OS": {"Homosapiens"}}], True) is None
