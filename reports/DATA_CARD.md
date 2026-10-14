# DATA_CARD — SV pathogenicity ablation dataset v1

**Intended use:** research comparison of short-read-only vs short+long-read information for SV pathogenicity classification in a controlled simulation. **Not for clinical use.**

## Composition
* 4,030 planted SVs (2,015 pathogenic / 2,015 benign), 3,766 DEL and 264 DUP, GRCh38 autosomes, 50 bp-1 Mb, one replicate, all heterozygous.
* Per SV: truth coordinates, label, group id (50 % reciprocal-overlap component, the unit for splitting), simulated short (ART) and HiFi reads, Manta/Delly/Sniffles2/cuteSV calls, 245 columns in `data/processed/features_all.parquet` (large intermediate BAMs/VCFs live outside the repo in `$DNA_DATA`).

## Provenance and label policy
| Source | Use | Version/date |
|---|---|---|
| ClinVar `variant_summary.txt.gz` | pathogenic (P/LP) and benign (B/LB), >=1 review star, DEL/DUP ("copy number loss/gain" mapped to DEL/DUP) | snapshot downloaded 2026-10 |
| dbVar ClinVar-SV `nstd102.GRCh38` | coordinate-interval widths (CI) for quality filtering only | 2026-06 |
| gnomAD-SV v4.1 sites | benign: PASS, AF >= 0.01, DEL/DUP | v4.1 |
| GENCODE v50, UCSC rmsk/segdup/simpleRepeat/gap/centromere, gnomAD v4.1 constraint | annotation features (no ClinGen/ClinVar-derived gene lists) | GRCh38 |
| Reference | Broad hg38 primary contigs chr1-22,X,Y,M | |

Excluded as features: allele frequency, source, review status, clinical text, cluster size, batch/replicate ids, coverage, truth-derived QC fields (enforced by `src/features/registry.py`).

## Processing
Filters and counts: see `truthset_manifest.json`. Conflicts (opposite labels, RO >= 0.5) removed; near-duplicates (RO >= 0.8) collapsed; at most 12 SVs per overlapping +-100 kb window; benign down-sampled to match pathogenic type x 0.25-decade size bins on common support, then both classes subsampled to about 2,000 each. Seeds: `configs/truthset.yaml`, per-SV seeds from `stable_seed(sv_id, ...)`.

## Known biases and limits
* Benign is 84 % gnomAD, pathogenic 100 % ClinVar; ClinVar pathogenic SVs are clinically ascertained (gene-rich), so gene overlap is nearly determinative (AUPRC ~0.98).
* Large benign CNVs are rare: few SVs >100 kb in the matched set.
* ClinVar array coordinates are imprecise; they are planted exactly as given.
* Simulated reads are cleaner than real reads (see `reports/data_audit/p9_realdata.md`).
* Short reads miss benign SVs at segmental-duplication breakpoints (selection in the both-detected cohort).

## Audits (all in `reports/data_audit/`)
`truthset_gate.md`, `sim_qc.md`, `matching.md`, `leakage.md`, `p9_realdata.md`.
