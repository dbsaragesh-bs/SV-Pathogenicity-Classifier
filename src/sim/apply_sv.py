"""Deterministic haplotype construction from a list of non-overlapping edits.

All coordinates are 0-based half-open in the frame of the input (reference) window.
An edit replaces seq[start:end] with `alt`.  Edits must not overlap; pure insertions
have start == end.  This one primitive expresses DEL, DUP (tandem), INV, INS and
small variants, so a single well-tested function builds every haplotype.
"""
from __future__ import annotations

from dataclasses import dataclass

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Edit:
    start: int
    end: int
    alt: str
    kind: str = "small"  # "sv" or "small"


def sv_edit(svtype: str, start: int, end: int, ref: str, ins_seq: str = "") -> Edit:
    """Edit for a structural variant on window coordinates [start, end)."""
    if svtype == "DEL":
        return Edit(start, end, "", "sv")
    if svtype == "DUP":  # tandem: second copy placed right after the first
        return Edit(end, end, ref[start:end], "sv")
    if svtype == "INV":
        return Edit(start, end, revcomp(ref[start:end]), "sv")
    if svtype == "INS":
        return Edit(start, start, ins_seq, "sv")
    raise ValueError(f"unsupported svtype {svtype}")


def apply_edits(seq: str, edits: list[Edit]) -> str:
    """Apply non-overlapping edits. Raises on overlap/out-of-range."""
    es = sorted(edits, key=lambda e: (e.start, e.end))
    prev_end = 0
    out = []
    pos = 0
    for e in es:
        if e.start < 0 or e.end > len(seq) or e.start > e.end:
            raise ValueError(f"edit out of range: {e}")
        if e.start < prev_end:
            raise ValueError(f"overlapping edits at {e}")
        out.append(seq[pos:e.start])
        out.append(e.alt)
        pos = e.end
        prev_end = e.end if e.end > e.start else e.start
        # an insertion at p followed by an edit starting at p is allowed
    out.append(seq[pos:])
    return "".join(out)
