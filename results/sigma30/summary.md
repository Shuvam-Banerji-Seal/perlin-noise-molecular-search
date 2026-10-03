# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=30.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## decane  (7 torsions, 32 atoms)

Reference: 1055 unique minima, 193 within window, E* = -6.299 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.351 ± 0.045 | 0.000 | 67.7 | 0.90 | 143 | 0.85 |
| perlin(L=4) | 0.319 ± 0.031 | 0.000 | 61.5 | 0.70 | 206 | 0.84 |
| gauss-corr(L=4) | 0.315 ± 0.022 | 0.823 | 60.8 | 0.40 | 162 | 0.85 |
| perlin(L=8) | 0.313 ± 0.034 | 0.000 | 60.5 | 0.80 | 103 | 0.86 |
| gauss-corr(L=8) | 0.294 ± 0.044 | 0.000 | 56.7 | 0.90 | 87 | 0.87 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.004 (Wilcoxon p = 0.898, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.020 (Wilcoxon p = 0.301, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.032 (Wilcoxon p = 0.084, n = 10)

## tetraglyme  (12 torsions, 37 atoms)

Reference: 3842 unique minima, 52 within window, E* = 88.954 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.454 ± 0.042 | 0.000 | 23.6 | 1.00 | 114 | 0.64 |
| gauss-corr(L=4) | 0.350 ± 0.039 | 0.000 | 18.2 | 1.00 | 139 | 0.70 |
| perlin(L=4) | 0.312 ± 0.068 | 0.000 | 16.2 | 0.90 | 112 | 0.71 |
| perlin(L=8) | 0.287 ± 0.055 | 0.000 | 14.9 | 1.00 | 124 | 0.68 |
| gauss-corr(L=8) | 0.277 ± 0.057 | 0.000 | 14.4 | 1.00 | 68 | 0.68 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.038 (Wilcoxon p = 0.141, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.010 (Wilcoxon p = 0.48, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.142 (Wilcoxon p = 0.00391, n = 10)

## ala3  (6 torsions, 42 atoms)

Reference: 388 unique minima, 8 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-iid | 0.487 ± 0.190 | 1.685 | 3.9 | 0.40 | 301 | 0.72 |
| gauss-corr(L=8) | 0.450 ± 0.251 | 0.000 | 3.6 | 0.60 | 205 | 0.64 |
| perlin(L=8) | 0.400 ± 0.211 | 2.643 | 3.2 | 0.20 | 280 | 0.66 |
| gauss-corr(L=4) | 0.388 ± 0.190 | 2.164 | 3.1 | 0.30 | 80 | 0.69 |
| perlin(L=4) | 0.375 ± 0.177 | 2.643 | 3.0 | 0.20 | 82 | 0.70 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.013 (Wilcoxon p = 1, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.050 (Wilcoxon p = 0.609, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.113 (Wilcoxon p = 0.438, n = 10)

## imatinib  (7 torsions, 68 atoms)

Reference: 1757 unique minima, 87 within window, E* = 64.812 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.236 ± 0.269 | 2.684 | 20.5 | 0.40 | 168 | 0.84 |
| perlin(L=4) | 0.224 ± 0.284 | 2.769 | 19.5 | 0.40 | 240 | 0.88 |
| perlin(L=8) | 0.154 ± 0.188 | 2.848 | 13.4 | 0.20 | 240 | 0.85 |
| gauss-corr(L=4) | 0.152 ± 0.234 | 2.926 | 13.2 | 0.30 | 57 | 0.86 |
| gauss-iid | 0.033 ± 0.083 | 4.951 | 2.9 | 0.10 | 488 | 0.92 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.072 (Wilcoxon p = 0.625, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.082 (Wilcoxon p = 0.285, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.191 (Wilcoxon p = 0.0234, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | +0.006 | 12 / 16 | 0.891 | 40 |
| gauss-corr(L=4) − gauss-iid | -0.030 | 11 / 21 | 0.0724 | 40 |
| perlin(L=4) − gauss-iid | -0.024 | 12 / 21 | 0.0983 | 40 |
| perlin(L=8) − gauss-corr(L=8) | -0.026 | 18 / 16 | 0.669 | 40 |
| gauss-corr(L=8) − gauss-iid | -0.017 | 14 / 22 | 0.307 | 40 |
| perlin(L=8) − gauss-iid | -0.043 | 11 / 23 | 0.0646 | 40 |
