"""Alignment-derived evidence features for a called SV (short and long reads).

Everything is computed from the BAM + reference only; no truth, no label, no simulator knobs.
Coordinates 0-based half-open on the *called* SV.
"""
from __future__ import annotations

import numpy as np
import pysam

NAN = float("nan")


def _depth_segments(bam, chrom, a, b, n_seg=16, seg_len=1000, min_pos=0):
    """Mean coverage over n_seg evenly spaced segments inside [a, b)."""
    L = b - a
    if L <= 0:
        return NAN
    if L <= seg_len * n_seg:
        segs = [(a, b)]
    else:
        step = L / n_seg
        segs = [(int(a + i * step), int(a + i * step) + seg_len) for i in range(n_seg)]
    tot = 0
    n = 0
    for s, e in segs:
        s = max(s, min_pos)
        if e <= s:
            continue
        cov = bam.count_coverage(chrom, s, e, quality_threshold=0, read_callback="all")
        tot += sum(int(x.sum()) for x in (np.asarray(c) for c in cov))
        n += e - s
    return tot / n if n else NAN


def depth_features(bam, chrom, start, end, clen) -> dict:
    L = end - start
    flank = int(min(max(L, 5000), 20000))
    gap = 300
    # interior: avoid the breakpoint zones
    ia, ib = (start + gap, end - gap) if L > 4 * gap else (start, end)
    d_in = _depth_segments(bam, chrom, ia, ib)
    d_l = _depth_segments(bam, chrom, max(start - gap - flank, 0), max(start - gap, 0), n_seg=8)
    d_r = _depth_segments(bam, chrom, min(end + gap, clen), min(end + gap + flank, clen), n_seg=8)
    fl = np.nanmean([d_l, d_r])
    ratio = d_in / fl if fl and fl > 0 else NAN
    return {
        "depth_in": d_in, "depth_flank": fl, "depth_ratio": ratio,
        "depth_log2ratio": float(np.log2(ratio)) if ratio and ratio > 0 else NAN,
        "depth_flank_asym": abs(d_l - d_r) / fl if fl and fl > 0 else NAN,
    }


# ------------------------------------------------------------------ short reads
def insert_stats(bam, chrom, a, b, n=20000):
    isz = []
    for r in bam.fetch(chrom, a, b):
        if r.is_proper_pair and r.is_read1 and not r.is_secondary and r.mapq >= 30:
            isz.append(abs(r.template_length))
            if len(isz) >= n:
                break
    if not isz:
        return 400.0, 50.0
    isz = np.array(isz)
    med = float(np.median(isz))
    mad = float(np.median(np.abs(isz - med)) * 1.4826)
    return med, max(mad, 10.0)


def short_features(bam: pysam.AlignmentFile, chrom, start, end, ins_mu, ins_sd, clen, w=500) -> dict:
    f = {}
    wins = {"l": (max(start - w, 0), start + w), "r": (max(end - w, 0), min(end + w, clen))}
    thr = ins_mu + 4 * ins_sd
    tot_reads = 0
    for side, (a, b) in wins.items():
        n = clip = split = disc = span = span_ok = 0
        q = []
        mapq0 = 0
        clip_pos = []
        other = wins["r" if side == "l" else "l"]
        for r in bam.fetch(chrom, a, b):
            if r.is_unmapped or r.is_secondary or r.is_duplicate:
                continue
            n += 1
            q.append(r.mapping_quality)
            mapq0 += r.mapping_quality < 10
            ct = r.cigartuples or []
            if ct:
                if ct[0][0] in (4, 5) and ct[0][1] >= 10 and a <= r.reference_start <= b:
                    clip += 1
                    clip_pos.append(r.reference_start)
                if ct[-1][0] in (4, 5) and ct[-1][1] >= 10 and a <= r.reference_end <= b:
                    clip += 1
                    clip_pos.append(r.reference_end)
            if r.has_tag("SA"):
                sa = r.get_tag("SA").split(";")[0].split(",")
                if sa[0] == chrom and other[0] - 200 <= int(sa[1]) <= other[1] + 200:
                    split += 1
            if r.is_paired and r.mate_is_mapped and not r.is_proper_pair:
                disc += 1
            if r.is_paired and r.mate_is_mapped and r.next_reference_name == chrom and r.is_read1 is not None:
                mp = r.next_reference_start
                if other[0] <= mp <= other[1] and side == "l":
                    span += 1
                    # DEL: left read forward, mate reverse;  DUP: left read reverse, mate forward
                    span_ok += int((not r.is_reverse) == r.mate_is_reverse)
        tot_reads += n
        f[f"s_{side}_reads"] = n
        f[f"s_{side}_clip"] = clip
        f[f"s_{side}_split_to_other"] = split
        f[f"s_{side}_disc"] = disc
        f[f"s_{side}_mapq_mean"] = float(np.mean(q)) if q else NAN
        f[f"s_{side}_lowmapq_frac"] = mapq0 / n if n else NAN
        f[f"s_{side}_clip_disp"] = float(np.std(clip_pos)) if len(clip_pos) > 1 else NAN
        if side == "l":
            f["s_span_pairs"] = span
            f["s_span_fr"] = span_ok
    f["s_clip_total"] = f["s_l_clip"] + f["s_r_clip"]
    f["s_split_total"] = f["s_l_split_to_other"] + f["s_r_split_to_other"]
    f["s_disc_total"] = f["s_l_disc"] + f["s_r_disc"]
    f["s_support_total"] = f["s_span_pairs"] + f["s_split_total"]
    f["s_support_frac"] = f["s_support_total"] / tot_reads if tot_reads else NAN
    f["s_mapq_mean"] = float(np.nanmean([f["s_l_mapq_mean"], f["s_r_mapq_mean"]]))
    f["s_lowmapq_frac"] = float(np.nanmean([f["s_l_lowmapq_frac"], f["s_r_lowmapq_frac"]]))
    f.update({"s_" + k: v for k, v in depth_features(bam, chrom, start, end, clen).items()})
    return f


# ------------------------------------------------------------------ long reads
def long_features(bam: pysam.AlignmentFile, chrom, start, end, clen, w=2000) -> dict:
    f = {}
    L = end - start
    pos_l: dict = {}
    pos_r: dict = {}
    meta: dict = {}
    ref_span = {"l": 0, "r": 0}
    n_tot = 0
    lowmq = 0
    lens = []
    nm_s = []
    tol = max(0.2 * L, 200)
    for side, p in (("l", start), ("r", end)):
        a, b = max(p - w, 0), min(p + w, clen)
        for r in bam.fetch(chrom, a, b):
            if r.is_unmapped or r.is_secondary:
                continue
            n_tot += 1
            lens.append(r.query_alignment_length)
            lowmq += r.mapping_quality < 20
            ct = r.cigartuples or []
            name = r.query_name
            supported = False
            # (1) clipped / split alignment ending at this breakpoint, joined to a supplementary alignment
            if r.has_tag("SA") and ct:
                if ct[-1][0] in (4, 5) and ct[-1][1] >= 200 and abs(r.reference_end - p) < w:
                    (pos_l if side == "l" else pos_r).setdefault(name, r.reference_end)
                    supported = True
                elif ct[0][0] in (4, 5) and ct[0][1] >= 200 and abs(r.reference_start - p) < w:
                    (pos_l if side == "l" else pos_r).setdefault(name, r.reference_start)
                    supported = True
            # (2) CIGAR deletion matching the SV interval (small/medium deletions)
            if not supported and L < 100_000:
                rp = r.reference_start
                for op, ln in ct:
                    if op in (0, 7, 8):
                        rp += ln
                    elif op == 2:
                        if ln >= 30 and abs(rp - start) < tol and abs(rp + ln - end) < tol:
                            pos_l.setdefault(name, rp)
                            pos_r.setdefault(name, rp + ln)
                            supported = True
                            break
                        rp += ln
                    elif op == 3:
                        rp += ln
            if supported:
                meta.setdefault(name, (r.mapping_quality, not r.is_reverse))
            elif r.reference_start < p - 100 and r.reference_end > p + 100:
                ref_span[side] += 1
            if r.has_tag("NM") and r.query_alignment_length:
                nm_s.append(r.get_tag("NM") / r.query_alignment_length)
    sup = set(pos_l) | set(pos_r)
    f["l_support_l"] = len(pos_l)
    f["l_support_r"] = len(pos_r)
    f["l_support_total"] = len(sup)
    f["l_ref_span_l"] = ref_span["l"]
    f["l_ref_span_r"] = ref_span["r"]
    den = len(sup) + (ref_span["l"] + ref_span["r"]) / 2
    f["l_allele_frac"] = len(sup) / den if den > 0 else NAN
    vl, vr = np.array(list(pos_l.values())), np.array(list(pos_r.values()))
    f["l_bp_disp_l"] = float(np.std(vl)) if len(vl) > 1 else NAN
    f["l_bp_disp_r"] = float(np.std(vr)) if len(vr) > 1 else NAN
    disp = [x for x in (f["l_bp_disp_l"], f["l_bp_disp_r"]) if x == x]
    f["l_bp_disp_max"] = float(max(disp)) if disp else NAN
    f["l_refined_start"] = float(np.median(vl)) if len(vl) else NAN
    f["l_refined_end"] = float(np.median(vr)) if len(vr) else NAN
    mq = [m[0] for m in meta.values()]
    st = [m[1] for m in meta.values()]
    f["l_sup_mapq_mean"] = float(np.mean(mq)) if mq else NAN
    f["l_sup_strand_bal"] = float(np.mean(st)) if st else NAN
    f["l_lowmapq_frac"] = lowmq / n_tot if n_tot else NAN
    f["l_readlen_mean"] = float(np.mean(lens)) if lens else NAN
    f["l_nm_rate"] = float(np.mean(nm_s)) if nm_s else NAN
    f["l_reads_at_bp"] = n_tot
    f.update({"l_" + k: v for k, v in depth_features(bam, chrom, start, end, clen).items()})
    return f


def microhomology(ref, chrom, start, end, maxshift=50) -> dict:
    """Breakpoint sequence homology for a deletion/duplication interval [start, end)."""
    try:
        c = ref[chrom]
        lo = max(start - maxshift - 1, 0)
        ls = str(c[lo:start]).upper()[::-1]            # bases going left from the start breakpoint
        le = str(c[end - (start - lo):end]).upper()[::-1]  # bases going left from the end breakpoint
        left = 0
        while left < min(len(ls), len(le), maxshift) and ls[left] == le[left]:
            left += 1
        rs = str(c[start:start + maxshift]).upper()
        re_ = str(c[end:end + maxshift]).upper()
        right = 0
        while right < min(len(rs), len(re_)) and rs[right] == re_[right]:
            right += 1
    except Exception:
        return {"mh_len": NAN, "mh_left": NAN, "mh_right": NAN}
    return {"mh_len": left + right, "mh_left": left, "mh_right": right}
