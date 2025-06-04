import sys
from pathlib import Path
import pytest

# Ensure project root is on sys.path so `translation` can be imported when
# running pytest from the tests directory.
sys.path.append(str(Path(__file__).resolve().parents[1]))

from translation import translate_function, biopython_translate_function

@pytest.fixture
def dna_examples():
    return {
        "seq1": ("ATGGCCTAA", "MA"),  # ATG->M, GCC->A, stop
        "seq2": ("ATTAGTTGA", "ISW"),  # ATT->I, AGT->S, TGA->W
        "seq3": ("ATAGGG", "MG"),  # ATA->M, GGG->G
        "seq4": ("AGGATG", ""),  # start with stop codon
    }


def test_translate_function(dna_examples):
    for seq, expected in dna_examples.values():
        assert translate_function(seq) == expected


def test_biopython_and_manual_same(dna_examples):
    for seq, _ in dna_examples.values():
        assert biopython_translate_function(seq) == translate_function(seq)
