# tests for motif loading and IUPAC → regex conversion
# test_motifs.py - Sanghyeon (Person 2)

import re

from src.motifs import load_motifs, iupac_to_regex


def test_iupac_dna_expansion_basic():
    # N should expand to any base
    regex = iupac_to_regex("ATNG", "dna")
    assert regex == "AT[ACGT]G"
    assert re.fullmatch(regex, "ATAG")
    assert re.fullmatch(regex, "ATCG")
    assert re.fullmatch(regex, "ATTG")
    assert re.fullmatch(regex, "ATGG")


def test_iupac_protein_expansion_basic():
    # X = any amino acid, Z = E or Q, B = D or N
    regex = iupac_to_regex("AXZB", "protein")
    # exact string is stable and easy to assert on
    assert regex == "A[ACDEFGHIKLMNPQRSTVWY][EQ][DN]"
    assert re.fullmatch(regex, "AWED")
    assert re.fullmatch(regex, "ANQD")
    assert re.fullmatch(regex, "AKEN")


def test_iupac_unknown_kind_raises():
    try:
        iupac_to_regex("ATN", "rna")
    except ValueError:
        # expected path
        pass
    else:
        assert False, "iupac_to_regex should raise ValueError for unknown kind"


def test_load_dna_motifs_from_json():
    motifs = load_motifs("dna")
    # JSON file defines three DNA motifs
    assert set(motifs.keys()) == {"TATA_box", "polyA_signal", "EcoRI_site"}

    tata = motifs["TATA_box"]
    assert hasattr(tata, "search")
    # simple sanity check: should match a classic TATA box
    seq = "GGGCTATAAATACCC"
    match = tata.search(seq)
    assert match is not None
    assert match.group(0).startswith("TATA")


def test_load_protein_motifs_from_json():
    motifs = load_motifs("protein")
    assert set(motifs.keys()) == {
        "N_glycosylation",
        "Proline_directed_phospho",
        "RGD_integrin_binding",
    }

    n_glyco = motifs["N_glycosylation"]
    # N-X-[ST]-X, where X != P
    assert n_glyco.search("NSTT") is not None
    assert n_glyco.search("NQSA") is not None
    # should NOT match if the second or fourth residue is proline
    assert n_glyco.search("NPST") is None
    assert n_glyco.search("NSTP") is None


def test_load_motifs_can_return_strings():
    motifs = load_motifs("dna", compiled=False)
    assert isinstance(motifs["EcoRI_site"], str)
    assert motifs["EcoRI_site"] == "GAATTC"


def test_load_motifs_invalid_kind_raises():
    try:
        load_motifs("rna")
    except ValueError:
        pass
    else:
        assert False, "load_motifs should raise ValueError for unknown kind"
