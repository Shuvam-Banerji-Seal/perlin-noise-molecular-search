# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=120.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## eicosane  (17 torsions, 62 atoms)

Reference: 39139 unique minima, 360 within window, E* = -7.495 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.069 ± 0.023 | 0.862 | 24.9 | 0.10 | 168 | 0.53 |
| perlin(L=8) | 0.009 ± 0.006 | 2.267 | 3.2 | 0.00 | - | 0.12 |
| gauss-corr(L=8) | 0.006 ± 0.006 | 2.567 | 2.2 | 0.00 | - | 0.15 |
| perlin-drift(L=8) | 0.004 ± 0.005 | 2.754 | 1.6 | 0.00 | - | 0.12 |
| gauss-ar1(L=8) | 0.004 ± 0.004 | 2.725 | 1.4 | 0.00 | - | 0.13 |
| gauss-iid | 0.003 ± 0.002 | 2.772 | 1.0 | 0.00 | - | 0.10 |
| gauss-anneal(L=8) | 0.003 ± 0.003 | 2.718 | 1.0 | 0.00 | - | 0.11 |
| perlin-anneal(L=8) | 0.002 ± 0.003 | 3.029 | 0.8 | 0.00 | - | 0.12 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.003 (Wilcoxon p = 0.371, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 0.0117, n = 10)

## hexaglyme  (18 torsions, 51 atoms)

Reference: 39103 unique minima, 176 within window, E* = 132.136 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.084 ± 0.049 | 0.962 | 14.7 | 0.10 | 132 | 0.45 |
| gauss-corr(L=8) | 0.005 ± 0.005 | 2.901 | 0.8 | 0.00 | - | 0.09 |
| perlin-drift(L=8) | 0.005 ± 0.005 | 3.399 | 0.8 | 0.00 | - | 0.07 |
| gauss-ar1(L=8) | 0.003 ± 0.006 | 3.219 | 0.5 | 0.00 | - | 0.10 |
| gauss-iid | 0.002 ± 0.003 | 3.360 | 0.3 | 0.00 | - | 0.08 |
| perlin(L=8) | 0.002 ± 0.004 | 3.914 | 0.3 | 0.00 | - | 0.09 |
| perlin-anneal(L=8) | 0.002 ± 0.004 | 3.812 | 0.3 | 0.00 | - | 0.08 |
| gauss-anneal(L=8) | 0.000 ± 0.000 | 4.036 | 0.0 | 0.00 | - | 0.11 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.003 (Wilcoxon p = 0.25, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.000 (Wilcoxon p = 1, n = 10)

## ala6  (12 torsions, 72 atoms)

Reference: 32204 unique minima, 4 within window, E* = 3.182 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.075 ± 0.121 | 4.020 | 0.3 | 0.00 | - | 0.02 |
| reset | 0.050 ± 0.105 | 4.677 | 0.2 | 0.00 | - | 0.26 |
| gauss-ar1(L=8) | 0.050 ± 0.105 | 3.842 | 0.2 | 0.10 | 478 | 0.02 |
| perlin(L=8) | 0.025 ± 0.079 | 4.664 | 0.1 | 0.00 | - | 0.02 |
| perlin-anneal(L=8) | 0.025 ± 0.079 | 4.460 | 0.1 | 0.00 | - | 0.02 |
| gauss-anneal(L=8) | 0.025 ± 0.079 | 4.020 | 0.1 | 0.00 | - | 0.02 |
| gauss-iid | 0.000 ± 0.000 | 4.447 | 0.0 | 0.00 | - | 0.02 |
| perlin-drift(L=8) | 0.000 ± 0.000 | 5.632 | 0.0 | 0.00 | - | 0.03 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.050 (Wilcoxon p = 0.625, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.025 (Wilcoxon p = 1, n = 10)

## methyl-oleate  (16 torsions, 57 atoms)

Reference: 38975 unique minima, 202 within window, E* = -6.646 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.036 ± 0.021 | 2.096 | 7.2 | 0.00 | - | 0.55 |
| gauss-corr(L=8) | 0.021 ± 0.016 | 2.249 | 4.2 | 0.10 | 211 | 0.13 |
| perlin-drift(L=8) | 0.011 ± 0.010 | 1.799 | 2.3 | 0.00 | - | 0.11 |
| perlin(L=8) | 0.008 ± 0.008 | 2.676 | 1.7 | 0.00 | - | 0.13 |
| perlin-anneal(L=8) | 0.007 ± 0.004 | 2.567 | 1.4 | 0.00 | - | 0.12 |
| gauss-iid | 0.006 ± 0.005 | 2.660 | 1.2 | 0.00 | - | 0.11 |
| gauss-ar1(L=8) | 0.006 ± 0.008 | 3.011 | 1.2 | 0.00 | - | 0.13 |
| gauss-anneal(L=8) | 0.005 ± 0.005 | 2.603 | 1.0 | 0.00 | - | 0.14 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.012 (Wilcoxon p = 0.156, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.002 (Wilcoxon p = 0.781, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=8) − gauss-corr(L=8) | -0.016 | 10 / 15 | 0.128 | 40 |
| gauss-corr(L=8) − gauss-iid | +0.024 | 19 / 6 | 0.000658 | 40 |
| perlin(L=8) − gauss-iid | +0.008 | 14 / 8 | 0.027 | 40 |
| gauss-iid − reset | -0.057 | 0 / 31 | 1.17e-06 | 40 |
