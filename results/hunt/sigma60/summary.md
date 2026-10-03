# Benchmark summary

seeds=10, steps=500 minimisations/run, sigma=60.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## eicosane  (17 torsions, 62 atoms)

Reference: 38373 unique minima, 986 within window, E* = -7.545 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.024 ± 0.008 | 0.911 | 23.5 | 0.10 | 168 | 0.53 |
| perlin(L=8) | 0.019 ± 0.006 | 1.185 | 18.3 | 0.20 | 360 | 0.27 |
| gauss-corr(L=8) | 0.016 ± 0.010 | 1.811 | 15.6 | 0.10 | 267 | 0.28 |
| gauss-iid | 0.010 ± 0.006 | 1.716 | 10.2 | 0.10 | 421 | 0.19 |
| gauss-anneal(L=8) | 0.010 ± 0.002 | 1.742 | 10.1 | 0.00 | - | 0.22 |
| gauss-ar1(L=8) | 0.009 ± 0.005 | 1.739 | 9.2 | 0.00 | - | 0.25 |
| perlin-anneal(L=8) | 0.009 ± 0.005 | 1.735 | 8.8 | 0.00 | - | 0.21 |
| perlin-drift(L=8) | 0.005 ± 0.002 | 2.217 | 5.0 | 0.00 | - | 0.25 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.003 (Wilcoxon p = 0.445, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.008 (Wilcoxon p = 0.0469, n = 10)

## hexaglyme  (18 torsions, 51 atoms)

Reference: 38203 unique minima, 57 within window, E* = 130.466 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-ar1(L=8) | 0.042 ± 0.068 | 2.580 | 2.4 | 0.10 | 220 | 0.13 |
| perlin(L=8) | 0.035 ± 0.034 | 2.295 | 2.0 | 0.20 | 399 | 0.13 |
| reset | 0.019 ± 0.023 | 2.632 | 1.1 | 0.00 | - | 0.45 |
| gauss-corr(L=8) | 0.018 ± 0.026 | 3.023 | 1.0 | 0.10 | 358 | 0.14 |
| gauss-anneal(L=8) | 0.016 ± 0.021 | 3.118 | 0.9 | 0.00 | - | 0.14 |
| perlin-drift(L=8) | 0.007 ± 0.012 | 3.535 | 0.4 | 0.00 | - | 0.15 |
| gauss-iid | 0.004 ± 0.007 | 3.911 | 0.2 | 0.00 | - | 0.14 |
| perlin-anneal(L=8) | 0.002 ± 0.006 | 4.066 | 0.1 | 0.00 | - | 0.13 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.018 (Wilcoxon p = 0.133, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.032 (Wilcoxon p = 0.00391, n = 10)

## ala6  (12 torsions, 72 atoms)

Reference: 28175 unique minima, 17 within window, E* = 3.182 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.053 ± 0.076 | 3.399 | 0.9 | 0.10 | 459 | 0.06 |
| gauss-ar1(L=8) | 0.035 ± 0.074 | 4.151 | 0.6 | 0.00 | - | 0.06 |
| perlin-anneal(L=8) | 0.029 ± 0.075 | 4.151 | 0.5 | 0.00 | - | 0.06 |
| gauss-anneal(L=8) | 0.018 ± 0.028 | 4.381 | 0.3 | 0.00 | - | 0.07 |
| perlin-drift(L=8) | 0.018 ± 0.040 | 4.649 | 0.3 | 0.00 | - | 0.05 |
| reset | 0.012 ± 0.025 | 4.677 | 0.2 | 0.00 | - | 0.26 |
| perlin(L=8) | 0.012 ± 0.025 | 4.664 | 0.2 | 0.00 | - | 0.06 |
| gauss-iid | 0.006 ± 0.019 | 4.020 | 0.1 | 0.00 | - | 0.07 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = -0.041 (Wilcoxon p = 0.156, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.006 (Wilcoxon p = 1, n = 10)

## methyl-oleate  (16 torsions, 57 atoms)

Reference: 37955 unique minima, 530 within window, E* = -6.691 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| perlin(L=8) | 0.025 ± 0.014 | 1.443 | 13.1 | 0.10 | 88 | 0.28 |
| gauss-ar1(L=8) | 0.015 ± 0.016 | 1.756 | 7.8 | 0.00 | - | 0.29 |
| gauss-corr(L=8) | 0.013 ± 0.006 | 1.830 | 7.1 | 0.00 | - | 0.33 |
| reset | 0.012 ± 0.007 | 2.140 | 6.6 | 0.00 | - | 0.55 |
| perlin-anneal(L=8) | 0.010 ± 0.010 | 1.950 | 5.5 | 0.00 | - | 0.26 |
| gauss-anneal(L=8) | 0.010 ± 0.010 | 2.286 | 5.2 | 0.00 | - | 0.25 |
| gauss-iid | 0.009 ± 0.008 | 2.154 | 4.9 | 0.00 | - | 0.26 |
| perlin-drift(L=8) | 0.008 ± 0.008 | 2.398 | 4.4 | 0.00 | - | 0.29 |

- perlin(L=8) − gauss-corr(L=8): Δcoverage = +0.011 (Wilcoxon p = 0.0742, n = 10)
- perlin(L=8) − gauss-iid: Δcoverage = +0.015 (Wilcoxon p = 0.0547, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=8) − gauss-corr(L=8) | -0.002 | 21 / 11 | 0.41 | 40 |
| gauss-corr(L=8) − gauss-iid | +0.018 | 20 / 7 | 0.00139 | 40 |
| perlin(L=8) − gauss-iid | +0.015 | 26 / 5 | 0.000224 | 40 |
| gauss-iid − reset | -0.010 | 5 / 22 | 0.000667 | 40 |
