# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=60.0 deg, kT=3.0 kcal/mol, window=3.0 kcal/mol

## decane  (7 torsions, 32 atoms)

Reference: 1719 unique minima, 193 within window, E* = -6.299 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=4) | 0.540 ± 0.011 | 0.823 | 104.3 | 0.40 | 340 | 0.73 |
| gauss-corr(L=4) | 0.538 ± 0.022 | 0.412 | 103.9 | 0.50 | 151 | 0.74 |
| gauss-iid | 0.523 ± 0.027 | 0.412 | 101.0 | 0.50 | 219 | 0.72 |
| gauss-corr(L=8) | 0.498 ± 0.045 | 0.000 | 96.1 | 0.70 | 134 | 0.76 |
| perlin(L=8) | 0.494 ± 0.044 | 0.412 | 95.4 | 0.50 | 172 | 0.75 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.002 (Wilcoxon p = 0.91, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.004 (Wilcoxon p = 0.674, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.017 (Wilcoxon p = 0.109, n = 10)

## tetraglyme  (12 torsions, 37 atoms)

Reference: 20510 unique minima, 46 within window, E* = 88.954 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.093 ± 0.038 | 1.597 | 4.3 | 0.10 | 380 | 0.61 |
| perlin(L=8) | 0.093 ± 0.054 | 1.494 | 4.3 | 0.20 | 420 | 0.62 |
| perlin(L=4) | 0.061 ± 0.055 | 1.634 | 2.8 | 0.20 | 308 | 0.63 |
| gauss-corr(L=4) | 0.052 ± 0.036 | 1.577 | 2.4 | 0.10 | 130 | 0.64 |
| gauss-iid | 0.028 ± 0.025 | 1.690 | 1.3 | 0.00 | - | 0.66 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.009 (Wilcoxon p = 0.906, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.000 (Wilcoxon p = 0.969, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.033 (Wilcoxon p = 0.125, n = 10)

## ala3  (6 torsions, 42 atoms)

Reference: 1020 unique minima, 9 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.867 ± 0.047 | 0.000 | 7.8 | 1.00 | 104 | 0.56 |
| perlin(L=4) | 0.867 ± 0.047 | 0.000 | 7.8 | 1.00 | 110 | 0.56 |
| perlin(L=8) | 0.867 ± 0.047 | 0.000 | 7.8 | 1.00 | 36 | 0.58 |
| gauss-iid | 0.856 ± 0.054 | 0.000 | 7.7 | 1.00 | 106 | 0.57 |
| gauss-corr(L=8) | 0.833 ± 0.094 | 0.000 | 7.5 | 0.80 | 56 | 0.58 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.000 (Wilcoxon p = 1, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.033 (Wilcoxon p = 0.438, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.011 (Wilcoxon p = 1, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | +0.004 | 10 / 10 | 0.881 | 30 |
| gauss-corr(L=4) − gauss-iid | +0.017 | 13 / 6 | 0.0668 | 30 |
| perlin(L=4) − gauss-iid | +0.020 | 15 / 5 | 0.0453 | 30 |
| perlin(L=8) − gauss-corr(L=8) | +0.010 | 11 / 12 | 0.843 | 30 |
| gauss-corr(L=8) − gauss-iid | +0.006 | 16 / 9 | 0.36 | 30 |
| perlin(L=8) − gauss-iid | +0.016 | 15 / 10 | 0.269 | 30 |
