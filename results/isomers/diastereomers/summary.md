# Benchmark summary

seeds=10, steps=300 minimisations/run, sigma=60.0 deg, kT=1.0 kcal/mol, window=3.0 kcal/mol

## LL-dipeptide  (4 torsions, 32 atoms)

Reference: 100 unique minima, 15 within window, E* = -10.400 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.953 ± 0.032 | 0.000 | 14.3 | 1.00 | 26 | 0.36 |
| gauss-iid | 0.927 ± 0.066 | 0.000 | 13.9 | 1.00 | 37 | 0.38 |
| gauss-corr(L=4) | 0.893 ± 0.064 | 0.000 | 13.4 | 1.00 | 15 | 0.45 |
| perlin(L=4) | 0.820 ± 0.100 | 0.000 | 12.3 | 1.00 | 18 | 0.42 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.073 (Wilcoxon p = 0.125, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.107 (Wilcoxon p = 0.0312, n = 10)

## LD-dipeptide  (4 torsions, 32 atoms)

Reference: 98 unique minima, 14 within window, E* = -10.893 kcal/mol

| method | coverage | gap (median) | n_low | hit rate | median hit step | acc. |
|---|---|---|---|---|---|---|
| reset | 0.957 ± 0.050 | 0.000 | 13.4 | 1.00 | 17 | 0.39 |
| gauss-corr(L=4) | 0.936 ± 0.063 | 0.000 | 13.1 | 1.00 | 12 | 0.45 |
| gauss-iid | 0.914 ± 0.081 | 0.000 | 12.8 | 1.00 | 22 | 0.44 |
| perlin(L=4) | 0.900 ± 0.069 | 0.000 | 12.6 | 1.00 | 8 | 0.41 |

- perlin(L=4) − gauss-corr(L=4): Δcoverage = -0.036 (Wilcoxon p = 0.562, n = 10)
- perlin(L=4) − gauss-iid: Δcoverage = -0.014 (Wilcoxon p = 0.844, n = 10)

## Pooled over molecules (paired by molecule and seed)

| comparison | Δcoverage | wins / losses | Wilcoxon p | n |
|---|---|---|---|---|
| perlin(L=4) − gauss-corr(L=4) | -0.055 | 3 / 10 | 0.0534 | 20 |
| gauss-corr(L=4) − gauss-iid | -0.006 | 6 / 8 | 0.636 | 20 |
| perlin(L=4) − gauss-iid | -0.060 | 4 / 11 | 0.0463 | 20 |
| gauss-iid − reset | -0.035 | 4 / 9 | 0.0615 | 20 |
