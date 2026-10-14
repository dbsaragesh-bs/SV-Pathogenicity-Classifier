import numpy as np
from src.sim.make_window import build_window

CFG = {"flank": 2000, "snv_density": 0.002, "indel_frac": 0.1}


class FakeFa(dict):
    pass


def fake_ref(n=20000, seed=0):
    rng = np.random.default_rng(seed)
    s = "".join(rng.choice(list("ACGT"), n))
    return {"chr1": s}


def sv(t, s=8000, e=10000):
    return dict(sv_id="T1", chrom="chr1", start=s, end=e, svtype=t)


def test_deletion_length_and_het():
    ref = fake_ref()
    w = build_window(sv("DEL"), ref, CFG)
    L = 2000 + 2000 + 2000
    # hap1 differs from hap2 by the SV length, up to net indel changes
    assert abs((len(w.hap1) - len(w.hap2)) - 2000) < 40


def test_dup_longer():
    w = build_window(sv("DUP"), fake_ref(), CFG)
    assert len(w.hap2) - len(w.hap1) > 1900


def test_hom_both_carry():
    w = build_window(sv("DEL"), fake_ref(), CFG, zygosity="hom")
    assert abs(len(w.hap1) - len(w.hap2)) < 40 and len(w.hap1) < 5000


def test_deterministic():
    a = build_window(sv("DEL"), fake_ref(), CFG)
    b = build_window(sv("DEL"), fake_ref(), CFG)
    assert a.hap1 == b.hap1 and a.hap2 == b.hap2
