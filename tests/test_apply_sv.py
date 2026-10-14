import pytest
from src.sim.apply_sv import Edit, apply_edits, revcomp, sv_edit

REF = "AAAACCCCGGGGTTTT"  # 16 bp


def test_del():
    assert apply_edits(REF, [sv_edit("DEL", 4, 8, REF)]) == "AAAAGGGGTTTT"


def test_dup_tandem():
    assert apply_edits(REF, [sv_edit("DUP", 4, 8, REF)]) == "AAAACCCCCCCCGGGGTTTT"


def test_inv():
    assert apply_edits(REF, [sv_edit("INV", 4, 12, REF)]) == "AAAA" + revcomp("CCCCGGGG") + "TTTT"


def test_ins():
    assert apply_edits(REF, [sv_edit("INS", 8, 8, REF, "ACGT")]) == "AAAACCCCACGTGGGGTTTT"


def test_small_plus_sv():
    out = apply_edits(REF, [Edit(0, 1, "T"), sv_edit("DEL", 4, 8, REF), Edit(14, 14, "AA")])
    assert out == "TAAAGGGGTTAATT"


def test_overlap_raises():
    with pytest.raises(ValueError):
        apply_edits(REF, [Edit(2, 6, ""), Edit(5, 7, "")])


def test_len_arith():
    for t, d in [("DEL", -4), ("DUP", 4), ("INV", 0)]:
        assert len(apply_edits(REF, [sv_edit(t, 4, 8, REF)])) == len(REF) + d
