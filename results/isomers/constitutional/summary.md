# Benchmark summary

seeds=10, steps=300 minimisations/run, sigma=60.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## n-decane  (7 torsions, 32 atoms)

Reference: 1431 unique minima, 193 within window, E* = -6.299 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.443 ± 0.042 | 0.000 | 85.5 | 0.80 | 146 | 0.42 |
| perlin(L=4) | 0.440 ± 0.023 | 0.000 | 84.9 | 0.70 | 240 | 0.47 |
| gauss-corr(L=4) | 0.433 ± 0.039 | 0.000 | 83.6 | 0.60 | 110 | 0.51 |
| reset | 0.427 ± 0.030 | 0.000 | 82.4 | 0.70 | 133 | 0.54 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.007 (Wilcoxon p = 0.645, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.003 (Wilcoxon p = 0.672, n = 10)

## 2-methylnonane  (6 torsions, 32 atoms)

Reference: 977 unique minima, 179 within window, E* = -0.651 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=4) | 0.512 ± 0.040 | 0.000 | 91.7 | 1.00 | 74 | 0.50 |
| gauss-iid | 0.506 ± 0.026 | 0.000 | 90.6 | 0.80 | 60 | 0.46 |
| reset | 0.498 ± 0.028 | 0.000 | 89.1 | 1.00 | 96 | 0.54 |
| gauss-corr(L=4) | 0.494 ± 0.035 | 0.000 | 88.5 | 1.00 | 99 | 0.52 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.018 (Wilcoxon p = 0.406, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 0.768, n = 10)

## 2,7-dimethyloctane  (5 torsions, 32 atoms)

Reference: 198 unique minima, 39 within window, E* = 4.992 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.928 ± 0.038 | 0.000 | 36.2 | 1.00 | 19 | 0.50 |
| reset | 0.928 ± 0.054 | 0.000 | 36.2 | 1.00 | 9 | 0.53 |
| gauss-corr(L=4) | 0.918 ± 0.036 | 0.000 | 35.8 | 1.00 | 26 | 0.55 |
| perlin(L=4) | 0.905 ± 0.024 | 0.000 | 35.3 | 1.00 | 24 | 0.52 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.013 (Wilcoxon p = 0.576, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.023 (Wilcoxon p = 0.352, n = 10)

## 4-propylheptane  (6 torsions, 32 atoms)

Reference: 396 unique minima, 25 within window, E* = 1.789 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.900 ± 0.051 | 0.000 | 22.5 | 1.00 | 62 | 0.33 |
| perlin(L=4) | 0.876 ± 0.095 | 0.000 | 21.9 | 1.00 | 24 | 0.40 |
| reset | 0.852 ± 0.060 | 0.000 | 21.3 | 1.00 | 61 | 0.40 |
| gauss-corr(L=4) | 0.848 ± 0.086 | 0.000 | 21.2 | 1.00 | 46 | 0.41 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.028 (Wilcoxon p = 0.438, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.024 (Wilcoxon p = 0.781, n = 10)

## 3,3-diethylhexane  (5 torsions, 32 atoms)

Reference: 171 unique minima, 3 within window, E* = 15.472 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 1.000 ± 0.000 | 0.000 | 3.0 | 1.00 | 38 | 0.28 |
| gauss-corr(L=4) | 1.000 ± 0.000 | 0.000 | 3.0 | 1.00 | 26 | 0.41 |
| perlin(L=4) | 0.967 ± 0.105 | 0.000 | 2.9 | 1.00 | 30 | 0.42 |
| reset | 0.933 ± 0.141 | 0.000 | 2.8 | 1.00 | 16 | 0.27 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.033 (Wilcoxon p = 1, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.033 (Wilcoxon p = 1, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | +0.001 | 22 / 17 | 0.375 | 50 |
| gauss-corr(L=4) − gauss-iid | -0.017 | 15 / 23 | 0.0555 | 50 |
| perlin(L=4) − gauss-iid | -0.015 | 16 / 18 | 0.494 | 50 |
| gauss-iid − reset | +0.028 | 25 / 11 | 0.0109 | 50 |
