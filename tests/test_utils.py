from src import utils


def test_clean_and_standardize_strings():
    assert utils.clean("  hello  ") == "hello"
    assert utils.clean("   ") is None
    assert utils.standard(" A B ") == "AB"
    assert utils.standard_upper(" a b ") == "AB"


def test_equiv_ignores_case_and_spaces():
    assert utils.equiv("Homo sapiens", "  homoSapiens ")
    assert not utils.equiv("Homo sapiens", "Pan troglodytes")
    assert not utils.equiv(None, "Pan troglodytes")


def test_margin_tolerance_and_matching_masses():
    assert utils.margin_tolerance(1000.0, 0.1) == 0.1
    assert utils.margin_tolerance(1000.0, 20) == 0.02

    assert utils.matching_masses(1000.0, 1000.05, 0.1)
    assert not utils.matching_masses(1000.0, 1000.2, 0.1)

    assert utils.matching_masses(1000.0, 1000.01, 20)
    assert not utils.matching_masses(1000.0, 1000.05, 20)


def test_is_ptm_accepts_expected_formats():
    assert utils.is_PTM(None, {"H", "D", "C"})
    assert utils.is_PTM("1H", {"H", "D", "C"})
    assert utils.is_PTM("0H", {"H", "D", "C"})
    assert utils.is_PTM("2H1D", {"H", "D", "C"})
    assert not utils.is_PTM("1X", {"H", "D", "C"})
