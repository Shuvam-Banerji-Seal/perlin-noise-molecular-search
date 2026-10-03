# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=60.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## decane  (7 torsions, 32 atoms)

Reference: 1715 unique minima, 193 within window, E* = -6.299 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| etkdg | 0.716 ± 0.022 | 0.000 | 138.2 | 1.00 | 92 | - |
| gauss-corr(L=2) | 0.616 ± 0.025 | 0.000 | 118.8 | 0.80 | 304 | 0.46 |
| perlin(L=4) | 0.616 ± 0.021 | 0.000 | 118.8 | 1.00 | 266 | 0.46 |
| perlin(L=2) | 0.612 ± 0.028 | 0.000 | 118.2 | 0.70 | 146 | 0.44 |
| gauss-iid | 0.607 ± 0.041 | 0.000 | 117.1 | 0.90 | 147 | 0.42 |
| gauss-corr(L=4) | 0.597 ± 0.031 | 0.000 | 115.3 | 0.80 | 150 | 0.50 |
| reset | 0.575 ± 0.047 | 0.000 | 111.0 | 0.80 | 144 | 0.55 |
| perlin(L=8) | 0.555 ± 0.055 | 0.000 | 107.2 | 0.80 | 172 | 0.54 |
| gauss-corr(L=8) | 0.541 ± 0.032 | 0.000 | 104.4 | 0.90 | 80 | 0.55 |

- perlin(L=2) − gauss-corr(L=2): Δcoverage = -0.003 (Wilcoxon p = 0.76, n = 10)
- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.018 (Wilcoxon p = 0.32, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.015 (Wilcoxon p = 0.264, n = 10)
- perlin(L=2) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 0.902, n = 10)

## tetraglyme  (12 torsions, 37 atoms)

Reference: 28615 unique minima, 56 within window, E* = 88.954 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.205 ± 0.065 | 0.000 | 11.5 | 0.70 | 196 | 0.26 |
| perlin(L=8) | 0.188 ± 0.044 | 0.000 | 10.5 | 0.90 | 165 | 0.24 |
| reset | 0.155 ± 0.116 | 1.494 | 8.7 | 0.20 | 126 | 0.46 |
| perlin(L=4) | 0.146 ± 0.059 | 1.450 | 8.2 | 0.40 | 256 | 0.23 |
| gauss-corr(L=4) | 0.120 ± 0.100 | 1.545 | 6.7 | 0.30 | 347 | 0.24 |
| gauss-iid | 0.089 ± 0.039 | 1.545 | 5.0 | 0.10 | 250 | 0.24 |
| gauss-corr(L=2) | 0.089 ± 0.061 | 1.494 | 5.0 | 0.30 | 372 | 0.21 |
| perlin(L=2) | 0.089 ± 0.024 | 1.545 | 5.0 | 0.00 | - | 0.19 |
| etkdg | 0.000 ± 0.000 | 4.644 | 0.0 | 0.00 | - | - |

- perlin(L=2) − gauss-corr(L=2): Δcoverage = -0.000 (Wilcoxon p = 0.992, n = 10)
- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.027 (Wilcoxon p = 0.447, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.018 (Wilcoxon p = 0.633, n = 10)
- perlin(L=2) − gauss-iid: Δcoverage = -0.000 (Wilcoxon p = 0.984, n = 10)

## ala3  (6 torsions, 42 atoms)

Reference: 834 unique minima, 8 within window, E* = -6.459 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=2) | 0.963 ± 0.084 | 0.000 | 7.7 | 1.00 | 51 | 0.22 |
| gauss-iid | 0.925 ± 0.105 | 0.000 | 7.4 | 1.00 | 126 | 0.25 |
| gauss-corr(L=8) | 0.875 ± 0.144 | 0.000 | 7.0 | 0.70 | 79 | 0.26 |
| reset | 0.838 ± 0.177 | 0.000 | 6.7 | 0.90 | 137 | 0.32 |
| perlin(L=2) | 0.838 ± 0.177 | 0.000 | 6.7 | 0.90 | 43 | 0.21 |
| perlin(L=8) | 0.838 ± 0.177 | 0.000 | 6.7 | 0.80 | 44 | 0.23 |
| perlin(L=4) | 0.762 ± 0.239 | 1.685 | 6.1 | 0.40 | 214 | 0.24 |
| gauss-corr(L=4) | 0.713 ± 0.213 | 0.000 | 5.7 | 0.60 | 66 | 0.25 |
| etkdg | 0.000 ± 0.000 | 5.666 | 0.0 | 0.00 | - | - |

- perlin(L=2) − gauss-corr(L=2): Δcoverage = -0.125 (Wilcoxon p = 0.0625, n = 10)
- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.050 (Wilcoxon p = 0.746, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.037 (Wilcoxon p = 0.562, n = 10)
- perlin(L=2) − gauss-iid: Δcoverage = -0.087 (Wilcoxon p = 0.172, n = 10)

## imatinib  (7 torsions, 68 atoms)

Reference: 3856 unique minima, 88 within window, E* = 64.812 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.582 ± 0.389 | 0.000 | 51.2 | 0.70 | 61 | 0.48 |
| gauss-corr(L=2) | 0.573 ± 0.383 | 0.000 | 50.4 | 0.70 | 108 | 0.52 |
| gauss-iid | 0.481 ± 0.395 | 0.000 | 42.3 | 0.60 | 148 | 0.53 |
| perlin(L=8) | 0.476 ± 0.393 | 0.000 | 41.9 | 0.60 | 175 | 0.45 |
| perlin(L=2) | 0.444 ± 0.372 | 0.000 | 39.1 | 0.60 | 227 | 0.49 |
| perlin(L=4) | 0.440 ± 0.440 | 1.385 | 38.7 | 0.50 | 37 | 0.48 |
| gauss-corr(L=4) | 0.349 ± 0.351 | 1.463 | 30.7 | 0.50 | 272 | 0.49 |
| reset | 0.268 ± 0.365 | 2.848 | 23.6 | 0.40 | 153 | 0.55 |
| etkdg | 0.112 ± 0.022 | 0.606 | 9.9 | 0.20 | 332 | - |

- perlin(L=2) − gauss-corr(L=2): Δcoverage = -0.128 (Wilcoxon p = 0.156, n = 10)
- perlin(L=4) − gauss-corr(L=4): Δcoverage = +0.091 (Wilcoxon p = 0.367, n = 10)
- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.106 (Wilcoxon p = 0.25, n = 10)
- perlin(L=2) − gauss-iid: Δcoverage = -0.036 (Wilcoxon p = 0.945, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=2) − gauss-corr(L=2) | -0.064 | 10 / 20 | 0.0161 | 40 |
| gauss-corr(L=2) − gauss-iid | +0.035 | 22 / 10 | 0.147 | 40 |
| perlin(L=2) − gauss-iid | -0.030 | 15 / 17 | 0.568 | 40 |
| perlin(L=4) − gauss-corr(L=4) | +0.046 | 22 / 14 | 0.162 | 40 |
| gauss-corr(L=4) − gauss-iid | -0.081 | 13 / 19 | 0.0273 | 40 |
| perlin(L=4) − gauss-iid | -0.034 | 21 / 12 | 0.844 | 40 |
| perlin(L=8) − gauss-corr(L=8) | -0.037 | 16 / 18 | 0.374 | 40 |
| gauss-corr(L=8) − gauss-iid | +0.025 | 16 / 16 | 0.432 | 40 |
| perlin(L=8) − gauss-iid | -0.011 | 17 / 17 | 0.959 | 40 |
| gauss-iid − reset | +0.066 | 24 / 9 | 0.0569 | 40 |
| etkdg − gauss-iid | -0.318 | 14 / 26 | 0.00271 | 40 |
