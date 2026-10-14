"""Window + haplotype construction with label-blind background variation."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from src.sim.apply_sv import Edit, apply_edits, sv_edit
from src.utils.common import stable_seed

BASES = "ACGT"


@dataclass
class Window:
    sv_id: str
    chrom: str
    ws: int            # window start on chrom (0-based)
    we: int
    sv_start: int      # SV in chrom coordinates, 0-based half-open
    sv_end: int
    svtype: str
    hap1: str
    hap2: str
    small: list = field(default_factory=list)  # (pos_in_window, ref, alt, gt1, gt2)


def background_variants(rng, ref: str, lo: int, hi: int, density: float, indel_frac: float,
                       p_gt=(0.4, 0.4, 0.2), margin: int = 50):
    """Random het/hom SNVs and short indels in window coordinates.

    Returns list of (pos, ref_allele, alt_allele, g1, g2). Variants avoid +-margin of the SV
    breakpoints and are >=10 bp apart so that edits never overlap.
    """
    L = len(ref)
    n = rng.poisson(L * density)
    pos = np.unique(rng.integers(0, L, n))
    if len(pos) > 1:
        keep = np.concatenate([[True], np.diff(pos) >= 10])
        pos = pos[keep]
    out = []
    gts = [(0, 1), (1, 0), (1, 1)]
    for p in pos:
        p = int(p)
        if abs(p - lo) < margin or abs(p - hi) < margin or p >= L - 10:
            continue
        r = ref[p]
        if r == "N":
            continue
        g = gts[rng.choice(3, p=p_gt)]
        if rng.random() >= indel_frac:
            alt = rng.choice([b for b in BASES if b != r])
            out.append((p, r, str(alt), g[0], g[1]))
        else:
            ln = int(min(5, rng.geometric(0.5)))
            if rng.random() < 0.5:  # deletion of ln bases (anchor-free representation)
                if "N" in ref[p:p + ln]:
                    continue
                out.append((p, ref[p:p + ln], "", g[0], g[1]))
            else:  # insertion before p
                ins = "".join(rng.choice(list(BASES), ln))
                out.append((p, "", ins, g[0], g[1]))
    return out


def _small_edit(v):
    p, r, a, _, _ = v
    return Edit(p, p + len(r), a, "small")


def build_window(sv: dict, ref_fa, cfg: dict, zygosity: str = "het", ins_seq: str = "") -> Window:
    chrom, s, e, t = sv["chrom"], int(sv["start"]), int(sv["end"]), sv["svtype"]
    flank = cfg["flank"]
    ws, we = s - flank, e + flank
    ref = str(ref_fa[chrom][ws:we]).upper()
    lo, hi = s - ws, e - ws
    rng = np.random.default_rng(stable_seed("bg", sv["sv_id"]))
    small = background_variants(rng, ref, lo, hi, cfg["snv_density"], cfg["indel_frac"])
    e_sv = sv_edit(t, lo, hi, ref, ins_seq)

    def hap_edits(h: int, carry_sv: bool):
        ed = []
        for v in small:
            if v[3 + h] != 1:
                continue
            if carry_sv and t in ("DEL", "INV") and lo - 1 <= v[0] < hi + 1:
                continue  # variant falls inside a deleted/inverted segment on this haplotype
            ed.append(_small_edit(v))
        if carry_sv:
            ed.append(e_sv)
        return ed

    sv_on = {"het": (False, True), "hom": (True, True)}[zygosity]
    hap1 = apply_edits(ref, hap_edits(0, sv_on[0]))
    hap2 = apply_edits(ref, hap_edits(1, sv_on[1]))
    return Window(sv["sv_id"], chrom, ws, we, s, e, t, hap1, hap2, small)
