# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 300 relaxations. Times are right-censored at the budget (counted as 301).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched Gaussian (baseline) | perlin(L=4) − gauss-corr(L=4) | -0.055 [-0.097, -0.013] | 0.0534 (0.14) | -2.5 [-11.9, +6.4] | 0.793 (0.79) |
| correlation vs none | gauss-corr(L=4) − gauss-iid | -0.006 [-0.047, +0.036] | 0.636 (0.64) | -9.2 [-22.1, +4.6] | 0.117 (0.35) |
| Perlin vs i.i.d. | perlin(L=4) − gauss-iid | -0.060 [-0.111, -0.010] | 0.0463 (0.14) | -11.7 [-22.3, +1.8] | 0.00692 (0.028) |
| Perlin vs MCMM reset | perlin(L=4) − reset | -0.095 [-0.139, -0.054] | 0.00283 (0.011) | -11.1 [-26.9, +4.0] | 0.126 (0.35) |

## LL-dipeptide: 15 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 1.00 | 28 | 0.820 | 0.660 | 0.847 | 0.182 | 1.000 | 0.032 |
| gauss-corr(L=4) | 1.00 | 30 | 0.893 | 0.713 | 0.913 | 0.117 | 1.000 | 0.026 |
| reset | 1.00 | 37 | 0.953 | 0.767 | 0.967 | 0.068 | 1.000 | 0.031 |
| gauss-iid | 1.00 | 39 | 0.927 | 0.780 | 0.940 | 0.088 | 1.000 | 0.031 |

## LD-dipeptide: 14 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 1.00 | 11 | 0.900 | 0.750 | 0.957 | 0.057 | 1.000 | 0.013 |
| gauss-corr(L=4) | 1.00 | 14 | 0.936 | 0.793 | 0.979 | 0.039 | 1.000 | 0.012 |
| gauss-iid | 1.00 | 23 | 0.914 | 0.771 | 0.950 | 0.091 | 1.000 | 0.015 |
| reset | 1.00 | 24 | 0.957 | 0.729 | 0.993 | 0.022 | 1.000 | 0.013 |
