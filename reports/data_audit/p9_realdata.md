# Real-data realism check (P9): simulated vs GIAB HG002 deletions

Real: 104 HG002 v5.0q GRCh38 benchmark DEL, size 52-6640 bp; Illumina 2x250 novoalign and Revio HiFi streamed from GIAB, subsampled to ~30x, run through the SAME Delly/Sniffles2/cuteSV and feature code (Manta crashed on the real novoalign BAM -> Delly only for short reads).
Simulated comparison set: 1444 benign DEL <= 6640 bp.

## Detection (recall of planted/benchmark SV)

| tech | simulated | real |
|---|---|---|
| short, Delly only | 0.850 | 0.106 |
| short, all callers | 0.918 | n/a |
| HiFi | 0.958 | 0.644 |

## short-read evidence features: real vs simulated

Adversarial validation AUROC (real vs simulated, 23 features, n_real=11, n_sim=1326): **0.970** (0.5 = indistinguishable).

Most shifted features (KS):

| feature       |    ks |   wass_std |   real_median |   sim_median |
|:--------------|------:|-----------:|--------------:|-------------:|
| s_l_mapq_mean | 0.909 |      3.245 |        70     |       59.755 |
| s_mapq_mean   | 0.909 |      3.16  |        70     |       59.713 |
| s_r_mapq_mean | 0.909 |      2.743 |        70     |       59.721 |
| s_r_reads     | 0.783 |      2.913 |       140     |      210     |
| s_l_reads     | 0.778 |      3.426 |       113     |      211     |
| s_r_clip_disp | 0.52  |      1.45  |       155.449 |       26.448 |
| s_l_disc      | 0.516 |      3.337 |         1     |       11     |
| s_disc_total  | 0.516 |      3.228 |         2     |       23     |
| s_r_disc      | 0.515 |      2.058 |         1     |       11     |
| s_l_clip_disp | 0.512 |      1.357 |       136.121 |       26.28  |

## long-read evidence features: real vs simulated

Adversarial validation AUROC (real vs simulated, 20 features, n_real=67, n_sim=1384): **0.981** (0.5 = indistinguishable).

Most shifted features (KS):

| feature           |    ks |   wass_std |   real_median |   sim_median |
|:------------------|------:|-----------:|--------------:|-------------:|
| l_ref_span_l      | 0.654 |      2.191 |         8     |       18     |
| l_ref_span_r      | 0.631 |      2.142 |         8     |       18     |
| l_readlen_mean    | 0.617 |      1.799 |     17500.3   |    16148.4   |
| l_bp_disp_r       | 0.557 |      1.114 |         8.7   |        0     |
| l_bp_disp_max     | 0.556 |      1.151 |        10.34  |        0     |
| l_bp_disp_l       | 0.555 |      1.956 |         4.789 |        0     |
| l_depth_log2ratio | 0.493 |      2.537 |        -1.452 |       -0.993 |
| l_depth_ratio     | 0.493 |      1.931 |         0.365 |        0.503 |
| l_allele_frac     | 0.491 |      2.325 |         0.667 |        0.5   |
| l_reads_at_bp     | 0.433 |      0.849 |        73     |       90     |

## Simulation-trained models applied to real HG002 benign-proxy deletions

- Arm A: 11 real SVs scored; predicted pathogenic (p>=0.5): 0.000; mean p 0.009
- Arm B: 11 real SVs scored; predicted pathogenic (p>=0.5): 0.000; mean p 0.005

Real HG002 SVs are germline variants of a healthy individual, so a low predicted-pathogenic fraction is expected. Real pathogenic SVs with reads are not available at scale; this check therefore tests false-positive behaviour and feature realism, not sensitivity.