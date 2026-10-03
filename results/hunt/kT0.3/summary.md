# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=60.0 deg, kT=0.3 kcal/mol, window=3.0 kcal/mol

## decane  (7 torsions, 32 atoms)

Reference: 1325 unique minima, 193 within window, E* = -6.299 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=4) | 0.647 ± 0.039 | 0.000 | 124.9 | 1.00 | 148 | 0.21 |
| gauss-iid | 0.615 ± 0.040 | 0.000 | 118.6 | 1.00 | 64 | 0.13 |
| gauss-corr(L=4) | 0.614 ± 0.031 | 0.000 | 118.5 | 1.00 | 70 | 0.23 |
| gauss-corr(L=8) | 0.477 ± 0.070 | 0.000 | 92.0 | 1.00 | 130 | 0.33 |
| perlin(L=8) | 0.460 ± 0.047 | 0.000 | 88.8 | 1.00 | 73 | 0.30 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.033 (Wilcoxon p = 0.0391, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.017 (Wilcoxon p = 0.711, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.033 (Wilcoxon p = 0.0625, n = 10)

## tetraglyme  (12 torsions, 37 atoms)

Reference: 12289 unique minima, 50 within window, E* = 88.954 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=4) | 0.256 ± 0.050 | 0.000 | 12.8 | 0.80 | 242 | 0.07 |
| gauss-corr(L=8) | 0.254 ± 0.043 | 0.000 | 12.7 | 1.00 | 210 | 0.14 |
| gauss-corr(L=4) | 0.236 ± 0.076 | 0.000 | 11.8 | 0.60 | 197 | 0.07 |
| perlin(L=8) | 0.212 ± 0.056 | 0.000 | 10.6 | 0.70 | 245 | 0.12 |
| gauss-iid | 0.206 ± 0.080 | 1.494 | 10.3 | 0.30 | 439 | 0.05 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.020 (Wilcoxon p = 0.57, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.042 (Wilcoxon p = 0.25, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = +0.050 (Wilcoxon p = 0.0977, n = 10)

## ala3  (6 torsions, 42 atoms)

Reference: 731 unique minima, 8 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=8) | 0.688 ± 0.159 | 0.843 | 5.5 | 0.50 | 8 | 0.19 |
| gauss-iid | 0.662 ± 0.177 | 1.685 | 5.3 | 0.40 | 124 | 0.18 |
| perlin(L=4) | 0.650 ± 0.165 | 0.843 | 5.2 | 0.50 | 31 | 0.17 |
| gauss-corr(L=8) | 0.613 ± 0.199 | 1.685 | 4.9 | 0.30 | 41 | 0.20 |
| gauss-corr(L=4) | 0.550 ± 0.105 | 1.685 | 4.4 | 0.40 | 46 | 0.19 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.100 (Wilcoxon p = 0.125, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.075 (Wilcoxon p = 0.312, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.013 (Wilcoxon p = 0.961, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | +0.051 | 18 / 7 | 0.0302 | 30 |
| gauss-corr(L=4) − gauss-iid | -0.028 | 12 / 13 | 0.382 | 30 |
| perlin(L=4) − gauss-iid | +0.023 | 17 / 11 | 0.368 | 30 |
| perlin(L=8) − gauss-corr(L=8) | +0.005 | 11 / 14 | 0.851 | 30 |
| gauss-corr(L=8) − gauss-iid | -0.047 | 8 / 16 | 0.0424 | 30 |
| perlin(L=8) − gauss-iid | -0.041 | 9 / 19 | 0.103 | 30 |
