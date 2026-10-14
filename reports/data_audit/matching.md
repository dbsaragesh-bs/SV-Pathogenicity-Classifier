# Truth matching and detection audit (P4)

Planted SVs: 4030 (2015 pathogenic)

## Detection rate by technology and matching stringency

| tech | loose RO>=0.3 | std RO>=0.5 | strict RO>=0.8 |
|---|---|---|---|
| short | 0.917 | 0.916 | 0.912 |
| hifi | 0.946 | 0.945 | 0.944 |

## Detection rate by label within size band (the detection-vs-label audit)

| size_band   |   ('detected_short', 0) |   ('detected_short', 1) |   ('detected_hifi', 0) |   ('detected_hifi', 1) |
|:------------|------------------------:|------------------------:|-----------------------:|-----------------------:|
| 50-300      |                   0.864 |                   0.957 |                  0.925 |                  0.935 |
| 300-1k      |                   0.975 |                   1     |                  0.951 |                  0.946 |
| 1k-5k       |                   0.989 |                   0.989 |                  0.96  |                  0.962 |
| 5k-30k      |                   0.684 |                   0.98  |                  0.937 |                  0.97  |
| >30k        |                   0.8   |                   0.983 |                  0.94  |                  0.958 |

## Detection rate by label within breakpoint repeat context

|                          |   ('detected_short', 'mean') |   ('detected_short', 'size') |   ('detected_hifi', 'mean') |   ('detected_hifi', 'size') |
|:-------------------------|-----------------------------:|-----------------------------:|----------------------------:|----------------------------:|
| ('SD breakpoint', 0)     |                        0.524 |                          374 |                       0.92  |                         374 |
| ('SD breakpoint', 1)     |                        0.763 |                           93 |                       0.978 |                          93 |
| ('unique breakpoint', 0) |                        0.933 |                         1641 |                       0.943 |                        1641 |
| ('unique breakpoint', 1) |                        0.986 |                         1922 |                       0.95  |                        1922 |

Interpretation: benign (gnomAD) SVs are enriched at segmental-duplication breakpoints (NAHR-mediated), where short-read callers lose recall. This is real detection biology, not a labelling error, but it makes C_both a selected cohort; analyses are repeated on C_all and on non-segdup strata.

## Breakpoint precision of matched calls (|called - true|, bp)

| tech | median start | median end | p90 | p99 | n>50bp |
|---|---|---|---|---|---|
| short | 1 | 0 | 2 | 50 | 72 |
| hifi | 1 | 0 | 2 | 9 | 31 |

## Calls not matching any planted SV (inside windows, DEL/DUP)

Total simulated window length 489 Mb.

| tech   |   false-positive calls |   per_Mb |
|:-------|-----------------------:|---------:|
| hifi   |                     55 |   0.1125 |
| short  |                     21 |   0.043  |