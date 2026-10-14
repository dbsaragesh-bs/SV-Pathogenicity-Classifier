"""P4 gate: detection stability, detection-vs-label audit, breakpoint precision, false positives."""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.utils.common import PROCESSED, REPORTS, ensure


def main():
    d = pd.read_parquet(PROCESSED / "features_all.parquet")
    fp = pd.read_parquet(PROCESSED / "features_all_fp.parquet")
    nz = pd.read_parquet(PROCESSED / "nuisance.parquet")
    L = ["# Truth matching and detection audit (P4)", "", f"Planted SVs: {len(d)} ({int(d.label.sum())} pathogenic)", "",
         "## Detection rate by technology and matching stringency", "", "| tech | loose RO>=0.3 | std RO>=0.5 | strict RO>=0.8 |", "|---|---|---|---|"]
    for t in ("short", "hifi"):
        L.append(f"| {t} | " + " | ".join(f"{d[f'det_{t}_{lv}'].mean():.3f}" for lv in ("loose", "std", "strict")) + " |")
    d["size_band"] = pd.cut(d.svlen, [0, 300, 1000, 5000, 30000, 1e7], labels=["50-300", "300-1k", "1k-5k", "5k-30k", ">30k"])
    d["segdup"] = np.where(d.o_ann_segdup_frac > 0.1, "segdup", "non-segdup")
    d["bp_segdup"] = np.where(d.o_ann_bp_segdup_max > 0.5, "SD breakpoint", "unique breakpoint")
    L += ["", "## Detection rate by label within size band (the detection-vs-label audit)", ""]
    g = d.groupby(["size_band", "label"], observed=True)[["detected_short", "detected_hifi"]].mean().round(3).unstack("label")
    L.append(g.to_markdown())
    L += ["", "## Detection rate by label within breakpoint repeat context", ""]
    g = d.groupby(["bp_segdup", "label"], observed=True)[["detected_short", "detected_hifi"]].agg(["mean", "size"]).round(3)
    L.append(g.to_markdown())
    L += ["", "Interpretation: benign (gnomAD) SVs are enriched at segmental-duplication breakpoints (NAHR-mediated), where "
          "short-read callers lose recall. This is real detection biology, not a labelling error, but it makes C_both a "
          "selected cohort; analyses are repeated on C_all and on non-segdup strata.", "",
          "## Breakpoint precision of matched calls (|called - true|, bp)", "", "| tech | median start | median end | p90 | p99 | n>50bp |", "|---|---|---|---|---|---|"]
    for t in ("short", "hifi"):
        e = pd.concat([d[f"{t}_bp_err_start"].abs(), d[f"{t}_bp_err_end"].abs()]).dropna()
        L.append(f"| {t} | {d[f'{t}_bp_err_start'].abs().median():.0f} | {d[f'{t}_bp_err_end'].abs().median():.0f} | {e.quantile(.9):.0f} | {e.quantile(.99):.0f} | {(e>50).sum()} |")
    mb = nz.win_len.sum() / 1e6
    L += ["", "## Calls not matching any planted SV (inside windows, DEL/DUP)", "",
          f"Total simulated window length {mb:.0f} Mb.", "",
          fp.groupby("tech").size().rename("false-positive calls").to_frame().assign(per_Mb=lambda x: (x.iloc[:, 0] / mb).round(4)).to_markdown()]
    ensure(REPORTS / "data_audit")
    (REPORTS / "data_audit" / "matching.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
