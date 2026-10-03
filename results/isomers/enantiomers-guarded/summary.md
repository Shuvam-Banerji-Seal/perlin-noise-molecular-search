# Benchmark summary

seeds=10, steps=300 minimisations/run, sigma=60.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## ala1-S  (2 torsions, 22 atoms)

Reference: 12 unique minima, 7 within window, E* = -17.732 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 1.000 ± 0.000 | 0.000 | 7.0 | 1.00 | 12 | 0.46 |
| gauss-iid | 0.971 ± 0.060 | 0.000 | 6.8 | 1.00 | 6 | 0.66 |
| gauss-corr(L=4) | 0.886 ± 0.113 | 0.000 | 6.2 | 1.00 | 4 | 0.66 |
| perlin(L=4) | 0.886 ± 0.060 | 0.000 | 6.2 | 1.00 | 4 | 0.64 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.000 (Wilcoxon p = 0.625, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.086 (Wilcoxon p = 0.0312, n = 10)

## ala1-R  (2 torsions, 22 atoms)

Reference: 12 unique minima, 7 within window, E* = -17.732 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 1.000 ± 0.000 | 0.000 | 7.0 | 1.00 | 12 | 0.45 |
| gauss-iid | 0.986 ± 0.045 | 0.000 | 6.9 | 1.00 | 3 | 0.66 |
| perlin(L=4) | 0.900 ± 0.069 | 0.000 | 6.3 | 1.00 | 4 | 0.64 |
| gauss-corr(L=4) | 0.814 ± 0.069 | 0.000 | 5.7 | 1.00 | 6 | 0.65 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.086 (Wilcoxon p = 0.0625, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.086 (Wilcoxon p = 0.0703, n = 10)

## ibuprofen-R  (4 torsions, 33 atoms)

Reference: 20 unique minima, 18 within window, E* = 23.790 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 1.000 ± 0.000 | 0.000 | 18.0 | 1.00 | 3 | 0.86 |
| gauss-corr(L=4) | 1.000 ± 0.000 | 0.000 | 18.0 | 1.00 | 3 | 0.88 |
| perlin(L=4) | 1.000 ± 0.000 | 0.000 | 18.0 | 1.00 | 2 | 0.88 |
| gauss-iid | 0.994 ± 0.018 | 0.000 | 17.9 | 1.00 | 2 | 0.88 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.000 (Wilcoxon p = n/a, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 1, n = 10)

## ibuprofen-S  (4 torsions, 33 atoms)

Reference: 19 unique minima, 18 within window, E* = 23.790 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 1.000 ± 0.000 | 0.000 | 18.0 | 1.00 | 1 | 0.86 |
| gauss-corr(L=4) | 1.000 ± 0.000 | 0.000 | 18.0 | 1.00 | 2 | 0.88 |
| perlin(L=4) | 1.000 ± 0.000 | 0.000 | 18.0 | 1.00 | 0 | 0.88 |
| gauss-iid | 0.994 ± 0.018 | 0.000 | 17.9 | 1.00 | 1 | 0.89 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.000 (Wilcoxon p = n/a, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 1, n = 10)

## propranolol-R  (6 torsions, 40 atoms)

Reference: 769 unique minima, 39 within window, E* = 72.925 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.662 ± 0.141 | 0.000 | 25.8 | 0.60 | 46 | 0.45 |
| reset | 0.582 ± 0.193 | 0.000 | 22.7 | 0.90 | 91 | 0.49 |
| perlin(L=4) | 0.567 ± 0.163 | 0.000 | 22.1 | 0.80 | 76 | 0.48 |
| gauss-corr(L=4) | 0.546 ± 0.197 | 0.000 | 21.3 | 0.80 | 63 | 0.46 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.021 (Wilcoxon p = 0.902, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.095 (Wilcoxon p = 0.467, n = 10)

## propranolol-S  (6 torsions, 40 atoms)

Reference: 771 unique minima, 39 within window, E* = 72.925 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.610 ± 0.187 | 0.000 | 23.8 | 0.70 | 30 | 0.50 |
| gauss-corr(L=4) | 0.587 ± 0.183 | 0.000 | 22.9 | 0.90 | 68 | 0.49 |
| gauss-iid | 0.579 ± 0.181 | 0.000 | 22.6 | 1.00 | 30 | 0.41 |
| perlin(L=4) | 0.479 ± 0.168 | 0.000 | 18.7 | 0.60 | 20 | 0.43 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.108 (Wilcoxon p = 0.184, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.100 (Wilcoxon p = 0.246, n = 10)

## ala3-SSS  (6 torsions, 42 atoms)

Reference: 728 unique minima, 8 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.775 ± 0.262 | 0.000 | 6.2 | 0.80 | 58 | 0.24 |
| reset | 0.738 ± 0.181 | 0.000 | 5.9 | 0.80 | 135 | 0.32 |
| gauss-corr(L=4) | 0.662 ± 0.221 | 0.000 | 5.3 | 0.60 | 66 | 0.26 |
| perlin(L=4) | 0.600 ± 0.211 | 1.685 | 4.8 | 0.20 | 56 | 0.25 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.062 (Wilcoxon p = 0.672, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.175 (Wilcoxon p = 0.156, n = 10)

## ala3-RRR  (6 torsions, 42 atoms)

Reference: 698 unique minima, 8 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.750 ± 0.186 | 0.000 | 6.0 | 0.90 | 36 | 0.22 |
| reset | 0.713 ± 0.177 | 0.000 | 5.7 | 0.90 | 61 | 0.31 |
| gauss-iid | 0.637 ± 0.181 | 0.000 | 5.1 | 0.80 | 29 | 0.24 |
| perlin(L=4) | 0.588 ± 0.156 | 1.685 | 4.7 | 0.40 | 35 | 0.25 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.163 (Wilcoxon p = 0.0312, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.050 (Wilcoxon p = 0.688, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | -0.028 | 19 / 27 | 0.155 | 80 |
| gauss-corr(L=4) − gauss-iid | -0.044 | 19 / 34 | 0.0292 | 80 |
| perlin(L=4) − gauss-iid | -0.073 | 17 / 37 | 0.00138 | 80 |
| gauss-iid − reset | -0.005 | 19 / 21 | 0.468 | 80 |
