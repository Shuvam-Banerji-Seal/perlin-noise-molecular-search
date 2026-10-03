# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=120.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## decane  (7 torsions, 32 atoms)

Reference: 1785 unique minima, 193 within window, E* = -6.299 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.507 ± 0.026 | 0.000 | 97.9 | 0.60 | 172 | 0.34 |
| perlin(L=4) | 0.493 ± 0.037 | 0.823 | 95.2 | 0.30 | 296 | 0.30 |
| gauss-corr(L=4) | 0.483 ± 0.046 | 0.823 | 93.2 | 0.20 | 336 | 0.30 |
| perlin(L=8) | 0.482 ± 0.028 | 0.412 | 93.0 | 0.50 | 237 | 0.33 |
| gauss-iid | 0.456 ± 0.031 | 0.823 | 88.0 | 0.30 | 262 | 0.29 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.010 (Wilcoxon p = 0.625, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.025 (Wilcoxon p = 0.0449, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.037 (Wilcoxon p = 0.123, n = 10)

## tetraglyme  (12 torsions, 37 atoms)

Reference: 22704 unique minima, 233 within window, E* = 90.403 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=8) | 0.039 ± 0.020 | 1.352 | 9.2 | 0.20 | 250 | 0.13 |
| gauss-corr(L=8) | 0.036 ± 0.022 | 1.392 | 8.5 | 0.00 | - | 0.16 |
| perlin(L=4) | 0.018 ± 0.006 | 1.694 | 4.1 | 0.00 | - | 0.14 |
| gauss-corr(L=4) | 0.017 ± 0.010 | 1.839 | 4.0 | 0.00 | - | 0.15 |
| gauss-iid | 0.012 ± 0.007 | 2.238 | 2.7 | 0.00 | - | 0.17 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.000 (Wilcoxon p = 0.879, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.003 (Wilcoxon p = 0.922, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 0.109, n = 10)

## ala3  (6 torsions, 42 atoms)

Reference: 844 unique minima, 8 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.900 ± 0.099 | 0.000 | 7.2 | 1.00 | 130 | 0.08 |
| gauss-corr(L=8) | 0.900 ± 0.079 | 0.000 | 7.2 | 1.00 | 146 | 0.11 |
| gauss-iid | 0.887 ± 0.109 | 0.000 | 7.1 | 1.00 | 73 | 0.08 |
| perlin(L=4) | 0.863 ± 0.092 | 0.000 | 6.9 | 1.00 | 89 | 0.07 |
| perlin(L=8) | 0.787 ± 0.103 | 0.000 | 6.3 | 0.60 | 31 | 0.07 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.037 (Wilcoxon p = 0.562, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.113 (Wilcoxon p = 0.0312, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.025 (Wilcoxon p = 0.75, n = 10)

## imatinib  (7 torsions, 68 atoms)

Reference: 3314 unique minima, 88 within window, E* = 64.812 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.533 ± 0.356 | 0.000 | 46.9 | 0.70 | 76 | 0.21 |
| perlin(L=4) | 0.532 ± 0.354 | 0.000 | 46.8 | 0.70 | 100 | 0.22 |
| perlin(L=8) | 0.434 ± 0.363 | 0.000 | 38.2 | 0.60 | 143 | 0.18 |
| gauss-corr(L=8) | 0.387 ± 0.386 | 1.463 | 34.1 | 0.50 | 85 | 0.20 |
| gauss-iid | 0.306 ± 0.367 | 2.926 | 26.9 | 0.40 | 124 | 0.21 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.001 (Wilcoxon p = 0.742, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.047 (Wilcoxon p = 0.523, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.226 (Wilcoxon p = 0.195, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | -0.007 | 16 / 19 | 0.471 | 40 |
| gauss-corr(L=4) − gauss-iid | +0.068 | 21 / 11 | 0.0189 | 40 |
| perlin(L=4) − gauss-iid | +0.061 | 20 / 11 | 0.108 | 40 |
| perlin(L=8) − gauss-corr(L=8) | -0.022 | 9 / 25 | 0.00866 | 40 |
| gauss-corr(L=8) − gauss-iid | +0.043 | 25 / 9 | 0.0204 | 40 |
| perlin(L=8) − gauss-iid | +0.021 | 23 / 13 | 0.245 | 40 |
