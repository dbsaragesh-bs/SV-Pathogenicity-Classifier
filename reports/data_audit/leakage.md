# Leakage and shortcut audit (P6)

All scores are out-of-fold over held-out chromosome groups. Chance AUROC = 0.5; chance AUPRC = prevalence.

| test | AUPRC | AUROC | prevalence | verdict |
|---|---|---|---|---|
| nuisance_only | 0.549 | 0.528 | 0.534 | PASS (≈chance) |
| size_only | 0.568 | 0.546 | 0.534 | PASS (≈chance) |
| type_only | 0.528 | 0.485 | 0.534 | PASS (≈chance) |
| chrom_only | 0.571 | 0.543 | 0.534 | PASS (≈chance) |
| position_only | 0.551 | 0.515 | 0.534 | PASS (≈chance) |
| n_callers_only | 0.548 | 0.530 | 0.534 | PASS (≈chance) |
| label_shuffle_armA | 0.541 | 0.513 | 0.534 | PASS (≈chance) |
| real_armA_reference | 0.986 | 0.987 | 0.534 | reference |
| adversarial_pipeline_cols | 0.562 | 0.530 | 0.534 | PASS (≈chance) |
| source_swap_train_af_test_clinvar | 0.978 | 0.890 | 0.877 | REVIEW |
| source_swap_train_clinvar_test_af | 0.977 | 0.964 | 0.577 | REVIEW |

Identical feature rows: 0; across different groups: 0

## Top Arm A features (sanity)

- a_ann_n_exons_hit: 0.331
- a_ann_max_cds_frac_affected: 0.187
- a_ann_cds_overlap_bp_log: 0.081
- a_ann_segdup_frac: 0.026
- a_ann_n_genes_disrupted: 0.024
- a_delly_present: 0.019
- a_ann_rep_ltr_frac: 0.016
- s_r_mapq_mean: 0.012
- a_svlen_log: 0.012
- s_mapq_mean: 0.011
- a_delly_qual: 0.011
- a_ann_rep_allrep_frac: 0.009
- a_ann_rep_l1_frac: 0.009
- a_ann_trf_frac: 0.009
- a_manta_sr_ref: 0.009