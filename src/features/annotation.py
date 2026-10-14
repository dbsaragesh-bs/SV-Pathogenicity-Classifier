"""Reference annotation tracks and per-interval annotation features.

All inputs are public, label-independent: GENCODE, gnomAD constraint (population-derived),
RepeatMasker, segmental duplications, simple repeats, assembly gaps/centromeres, reference GC.
Coordinates: 0-based half-open everywhere.
"""
from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd

from src.utils.common import RAW


class MergedTrack:
    """Union of intervals per chromosome with O(log n) overlap-length queries."""

    def __init__(self, df: pd.DataFrame):
        self.s, self.e, self.cum = {}, {}, {}
        for c, g in df.groupby("chrom"):
            iv = g[["start", "end"]].to_numpy(np.int64)
            iv = iv[np.argsort(iv[:, 0], kind="stable")]
            ms, me = [iv[0, 0]], [iv[0, 1]]
            for a, b in iv[1:]:
                if a <= me[-1]:
                    me[-1] = max(me[-1], b)
                else:
                    ms.append(a)
                    me.append(b)
            ms, me = np.array(ms), np.array(me)
            self.s[c], self.e[c] = ms, me
            self.cum[c] = np.concatenate([[0], np.cumsum(me - ms)])

    def _cover(self, c, x):
        """Total covered bp in [0, x)."""
        s, e, cum = self.s[c], self.e[c], self.cum[c]
        i = np.searchsorted(s, x, side="right")  # intervals with start <= x
        if i == 0:
            return 0
        base = cum[i - 1]
        return base + min(x, e[i - 1]) - s[i - 1]

    def overlap_bp(self, c, a, b) -> int:
        if c not in self.s or b <= a:
            return 0
        return int(self._cover(c, b) - self._cover(c, a))

    def frac(self, c, a, b) -> float:
        return self.overlap_bp(c, a, b) / max(b - a, 1)


def _read_ucsc(path, cols):
    return pd.read_csv(path, sep="\t", header=None, usecols=list(cols.keys()), names=None,
                       dtype=str, compression="infer").rename(columns=cols)


class Annotations:
    def __init__(self, ref_fa=None, load_gencode=True):
        self.ref = ref_fa
        self._load_repeats()
        self._load_gaps()
        if load_gencode:
            self._load_genes()

    # ------------------------------------------------------------------ loading
    def _load_repeats(self):
        rm = pd.read_csv(RAW / "rmsk.txt.gz", sep="\t", header=None, usecols=[5, 6, 7, 10, 11, 12],
                         names=["chrom", "start", "end", "name", "cls", "fam"])
        rm["start"] = rm.start.astype(int)
        rm["end"] = rm.end.astype(int)
        self.rep = {
            "alu": MergedTrack(rm[rm.fam == "Alu"]),
            "l1": MergedTrack(rm[rm.fam == "L1"]),
            "ltr": MergedTrack(rm[rm.cls == "LTR"]),
            "satellite": MergedTrack(rm[rm.cls == "Satellite"]),
            "simple": MergedTrack(rm[rm.cls.isin(["Simple_repeat", "Low_complexity"])]),
            "allrep": MergedTrack(rm),
        }
        sd = pd.read_csv(RAW / "genomicSuperDups.txt.gz", sep="\t", header=None, usecols=[1, 2, 3, 26],
                         names=["chrom", "start", "end", "fracMatch"])
        self.segdup = MergedTrack(sd)
        self.segdup_hi = MergedTrack(sd[sd.fracMatch >= 0.98])
        sr = pd.read_csv(RAW / "simpleRepeat.txt.gz", sep="\t", header=None, usecols=[1, 2, 3],
                         names=["chrom", "start", "end"])
        self.trf = MergedTrack(sr)

    def _load_gaps(self):
        g = pd.read_csv(RAW / "gap.txt.gz", sep="\t", header=None, usecols=[1, 2, 3, 7],
                        names=["chrom", "start", "end", "type"])
        c = pd.read_csv(RAW / "centromeres.txt.gz", sep="\t", header=None, usecols=[1, 2, 3],
                        names=["chrom", "start", "end"])
        both = pd.concat([g[["chrom", "start", "end"]], c], ignore_index=True)
        self.gap = MergedTrack(both)
        self.gap_s = {ch: gg.sort_values("start")[["start", "end"]].to_numpy() for ch, gg in both.groupby("chrom")}

    def _load_genes(self):
        genes, tx = {}, {}
        rows_exon, rows_cds = [], []
        with gzip.open(RAW / next(p.name for p in RAW.glob("gencode.v*.annotation.gtf.gz")), "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                ft = f[2]
                if ft not in ("gene", "transcript", "exon", "CDS"):
                    continue
                attr = f[8]
                gid = re.search(r'gene_id "([^"]+)"', attr).group(1)
                if ft == "gene":
                    gt = re.search(r'gene_type "([^"]+)"', attr).group(1)
                    gn = re.search(r'gene_name "([^"]+)"', attr).group(1)
                    genes[gid] = (f[0], int(f[3]) - 1, int(f[4]), f[6], gt, gn)
                elif ft == "transcript":
                    tid = re.search(r'transcript_id "([^"]+)"', attr).group(1)
                    pri = 2 if "MANE_Select" in attr else (1 if "Ensembl_canonical" in attr else 0)
                    tx[tid] = (gid, pri, int(f[4]) - int(f[3]) + 1)
                else:
                    tid = re.search(r'transcript_id "([^"]+)"', attr).group(1)
                    (rows_exon if ft == "exon" else rows_cds).append((tid, f[0], int(f[3]) - 1, int(f[4])))
        # one representative transcript per gene: MANE > canonical > longest
        best = {}
        for tid, (gid, pri, ln) in tx.items():
            k = (pri, ln)
            if gid not in best or k > best[gid][0]:
                best[gid] = (k, tid)
        rep_t = {v[1] for v in best.values()}
        gdf = pd.DataFrame([(g, *v) for g, v in genes.items()],
                           columns=["gid", "chrom", "start", "end", "strand", "gtype", "gname"])
        gdf = gdf[gdf.chrom.isin([f"chr{i}" for i in range(1, 23)])]
        gdf = gdf.merge(pd.DataFrame(self._constraint(), columns=["gname", "loeuf", "pli"]), on="gname", how="left")
        self.genes = {c: g.reset_index(drop=True) for c, g in gdf.groupby("chrom")}
        ex = pd.DataFrame([r for r in rows_exon if r[0] in rep_t], columns=["tid", "chrom", "start", "end"])
        cd = pd.DataFrame([r for r in rows_cds if r[0] in rep_t], columns=["tid", "chrom", "start", "end"])
        ex["gid"] = ex.tid.map(lambda t: tx[t][0])
        cd["gid"] = cd.tid.map(lambda t: tx[t][0])
        self.exons = {c: g.reset_index(drop=True) for c, g in ex.groupby("chrom")}
        self.cds = MergedTrack(cd)
        self.exon_track = MergedTrack(ex)
        pc = gdf[gdf.gtype == "protein_coding"]
        self.tss = {c: np.sort(np.where(g.strand == "+", g.start, g.end)) for c, g in pc.groupby("chrom")}
        self.cds_len = cd.assign(l=cd.end - cd.start).groupby("gid").l.sum().to_dict()
        self._cds_by_gene = {g: x[["start", "end"]].to_numpy() for g, x in cd.groupby("gid")}

    @staticmethod
    def _constraint():
        p = RAW / "gnomad.v4.1.constraint_metrics.tsv"
        d = pd.read_csv(p, sep="\t", low_memory=False)
        d = d[(d.get("canonical", True) == True) | (d.get("mane_select", False) == True)]  # noqa: E712
        col_l = "lof.oe_ci.upper"
        col_p = "lof.pLI"
        d = d[["gene", col_l, col_p]].dropna(subset=["gene"]).drop_duplicates("gene")
        return list(map(tuple, d.to_numpy()))

    # ------------------------------------------------------------------ features
    def gc(self, c, a, b):
        if self.ref is None or b <= a:
            return np.nan
        s = str(self.ref[c][max(a, 0):b]).upper()
        n = len(s) - s.count("N")
        return (s.count("G") + s.count("C")) / n if n else np.nan

    def interval_features(self, c: str, a: int, b: int, prefix: str = "") -> dict:
        """Annotation features for one interval [a, b)."""
        a, b = int(a), int(b)
        L = max(b - a, 1)
        f: dict = {}
        # repeats / segdups: SV body and +-1kb breakpoint flanks
        for name, tr in self.rep.items():
            f[f"rep_{name}_frac"] = tr.frac(c, a, b)
        f["trf_frac"] = self.trf.frac(c, a, b)
        f["segdup_frac"] = self.segdup.frac(c, a, b)
        f["segdup_hi_frac"] = self.segdup_hi.frac(c, a, b)
        for side, (x0, x1) in {"l": (a - 1000, a + 1000), "r": (b - 1000, b + 1000)}.items():
            f[f"bp_{side}_segdup"] = self.segdup.frac(c, x0, x1)
            f[f"bp_{side}_repeat"] = self.rep["allrep"].frac(c, x0, x1)
            f[f"bp_{side}_gc"] = self.gc(c, x0, x1)
        f["bp_segdup_both"] = float(f["bp_l_segdup"] > 0.5 and f["bp_r_segdup"] > 0.5)
        f["bp_repeat_max"] = max(f["bp_l_repeat"], f["bp_r_repeat"])
        f["bp_segdup_max"] = max(f["bp_l_segdup"], f["bp_r_segdup"])
        f["bp_gc_diff"] = abs(f["bp_l_gc"] - f["bp_r_gc"])
        f["gc_body"] = self.gc(c, a, b) if L <= 2_000_000 else np.nan
        f["gap_overlap_frac"] = self.gap.frac(c, a, b)
        gs = self.gap_s.get(c)
        if gs is not None and len(gs):
            d = np.minimum(np.abs(gs[:, 0] - b), np.abs(gs[:, 1] - a))
            inside = (gs[:, 0] < b) & (gs[:, 1] > a)
            f["dist_gap_log"] = float(np.log1p(0 if inside.any() else d.min()))
        else:
            f["dist_gap_log"] = np.nan
        clen = len(self.ref[c]) if self.ref is not None else np.nan
        f["dist_chrom_end_log"] = float(np.log1p(min(a, clen - b))) if clen == clen else np.nan

        # genes
        g = self.genes.get(c)
        ov = g[(g.start < b) & (g.end > a)] if g is not None else pd.DataFrame()
        pcg = ov[ov.gtype == "protein_coding"] if len(ov) else ov
        f["n_genes_pc"] = len(pcg)
        f["n_genes_nc"] = len(ov) - len(pcg)
        f["n_genes_contained"] = int(((pcg.start >= a) & (pcg.end <= b)).sum()) if len(pcg) else 0
        f["n_genes_disrupted"] = int(((pcg.start < a) & (pcg.end > a)).sum() + ((pcg.start < b) & (pcg.end > b)).sum()) if len(pcg) else 0
        f["cds_overlap_bp_log"] = float(np.log1p(self.cds.overlap_bp(c, a, b)))
        f["exon_overlap_bp_log"] = float(np.log1p(self.exon_track.overlap_bp(c, a, b)))
        exs = self.exons.get(c)
        if exs is not None and len(pcg):
            e = exs[exs.gid.isin(pcg.gid)]
            touched = e[(e.start < b) & (e.end > a)]
            f["n_exons_hit"] = len(touched)
            f["n_exons_full"] = int(((touched.start >= a) & (touched.end <= b)).sum())
            f["bp_in_exon"] = int(((e.start < a) & (e.end > a)).any()) + int(((e.start < b) & (e.end > b)).any())
            # max fraction of any gene's CDS affected
            fr = []
            for gid in pcg.gid:
                tot = self.cds_len.get(gid, 0)
                if tot:
                    sub = self.cds_df_gene(gid, c)
                    o = np.clip(np.minimum(sub[:, 1], b) - np.maximum(sub[:, 0], a), 0, None).sum()
                    fr.append(o / tot)
            f["max_cds_frac_affected"] = max(fr) if fr else 0.0
        else:
            f.update(n_exons_hit=0, n_exons_full=0, bp_in_exon=0, max_cds_frac_affected=0.0)
        if len(pcg):
            lo = pcg.loeuf.dropna()
            f["min_loeuf"] = float(lo.min()) if len(lo) else np.nan
            f["mean_loeuf"] = float(lo.mean()) if len(lo) else np.nan
            f["max_pli"] = float(pcg.pli.max()) if pcg.pli.notna().any() else np.nan
            f["n_constrained"] = int((lo < 0.6).sum())
            f["sum_pli"] = float(pcg.pli.fillna(0).sum())
        else:
            f.update(min_loeuf=np.nan, mean_loeuf=np.nan, max_pli=np.nan, n_constrained=0, sum_pli=0.0)
        ts = self.tss.get(c)
        if ts is not None and len(ts):
            i = np.searchsorted(ts, [a, b])
            cand = np.concatenate([ts[np.clip(i - 1, 0, len(ts) - 1)], ts[np.clip(i, 0, len(ts) - 1)]])
            inside = ((ts >= a) & (ts <= b)).any()
            f["dist_tss_log"] = float(np.log1p(0 if inside else np.min(np.minimum(np.abs(cand - a), np.abs(cand - b)))))
        else:
            f["dist_tss_log"] = np.nan
        f["svlen_log"] = float(np.log10(L))
        return {prefix + k: v for k, v in f.items()}

    def cds_df_gene(self, gid, c):
        return self._cds_by_gene[gid]
