# Simulation QC (P2 gate)

Windows simulated: 4030

## Achieved vs target coverage

| tech | median ratio achieved/target | p5 | p95 |
|---|---|---|---|
| short | 1.087 | 0.986 | 1.139 |
| hifi | 1.005 | 1.001 | 1.009 |

## Label-blindness: simulation parameters, pathogenic vs benign (KS test)

| parameter | median path | median benign | KS stat | p |
|---|---|---|---|---|
| cov_short_target | 30.12 | 30.12 | 0.004 | 1 |
| cov_long_target | 30.7 | 30.7 | 0.003 | 1 |
| n_small | 103 | 103 | 0.027 | 0.464 |
| short_cov_ach | 32.37 | 32.3 | 0.031 | 0.278 |
| hifi_cov_ach | 30.83 | 30.83 | 0.011 | 1 |

Parameters failing p<0.01: none
(`n_small` scales with window length = SV size + flank; size is matched between classes, so it should pass too.)