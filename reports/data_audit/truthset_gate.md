# Truth-set gate (P1)

Matched set: 4030 SVs (2015 pathogenic).

## Schema / integrity

- no_null_keys: True
- start_lt_end: True
- svlen_consistent: True
- chrom_allowed: True
- label_values: True
- unique_sv_id: True
- cross_label_conflicts: 0
- near_duplicates_ro80: 0

Size distribution KS (log10 length): D=0.028, p=0.418

## Trivial baselines (held-out chromosome groups, expect ≈ chance)

| baseline | AUPRC | AUROC |
|---|---|---|
| size_only | 0.504 | 0.512 |
| type_only | 0.489 | 0.484 |
| chrom_only | 0.523 | 0.527 |
| position_only | 0.510 | 0.501 |
| size+type | 0.507 | 0.514 |

## Type x label

| svtype   |    0 |    1 |
|:---------|-----:|-----:|
| DEL      | 1883 | 1883 |
| DUP      |  132 |  132 |

## Source x label

| source    |    0 |    1 |
|:----------|-----:|-----:|
| clinvar   |  317 | 2015 |
| gnomad_sv | 1698 |    0 |