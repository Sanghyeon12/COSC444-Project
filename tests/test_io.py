import io
import pytest

from src.io import read_fasta, parse_fasta, detect_alphabet, DNA, PROTEIN, validate_sequence


def test_parse_simple_fasta():
    fasta = """>seq1 desc
    ACGT
    ACGT
    >seq2
    MKWVTF
    """
    records = list(parse_fasta(io.StringIO(fasta)))
    assert len(records) == 2
    assert records[0].id == "seq1"
    assert records[0].description == "desc"
    assert records[0].sequence == "ACGTACGT"
    assert records[1].id == "seq2"
    assert records[1].sequence == "MKWVTF"


def test_read_fasta_from_handle():
    fasta = """>a
    ACGT
    """
    records = read_fasta(io.StringIO(fasta))
    assert records[0].id == "a"
    assert records[0].sequence == "ACGT"


def test_invalid_fasta_raises():
    fasta = """ACGT
    >seq
    ACGT
    """
    with pytest.raises(ValueError):
        list(parse_fasta(io.StringIO(fasta)))


def test_detect_alphabet():
    assert detect_alphabet("ACGTACGT") == "DNA"
    assert detect_alphabet("MKWVTFIS") == "PROTEIN"
    assert detect_alphabet("ACGTMK") == "UNKNOWN"


def test_validate_sequence():
    validate_sequence("ACGT", DNA)
    validate_sequence("MKWVTF", PROTEIN)

    with pytest.raises(ValueError):
        validate_sequence("ACGTB", DNA)

    with pytest.raises(ValueError):
        validate_sequence("MKWVTFZ", PROTEIN)
