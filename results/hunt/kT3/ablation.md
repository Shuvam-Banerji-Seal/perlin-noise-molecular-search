# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 500 relaxations. Times are right-censored at the budget (counted as 501).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched (L=4), kT=3 | perlin(L=4) − gauss-corr(L=4) | +0.003 [-0.018, +0.024] | 0.97 (1) | +31.7 [-23.8, +89.7] | 0.391 (1) |
| Perlin vs matched (L=8), kT=3 | perlin(L=8) − gauss-corr(L=8) | +0.007 [-0.019, +0.035] | 0.831 (1) | -3.9 [-81.3, +72.7] | 0.884 (1) |
| correlation vs none, kT=3 | gauss-corr(L=8) − gauss-iid | +0.001 [-0.036, +0.033] | 0.427 (1) | -3.7 [-79.2, +74.3] | 0.695 (1) |

## decane: 193 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 0.70 | 305 | 0.498 | 0.138 | 0.543 | 0.283 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.50 | 348 | 0.538 | 0.154 | 0.577 | 0.260 | 1.000 | 0.000 |
| perlin(L=8) | 0.50 | 356 | 0.494 | 0.122 | 0.537 | 0.285 | 1.000 | 0.000 |
| gauss-iid | 0.50 | 360 | 0.523 | 0.148 | 0.559 | 0.269 | 1.000 | 0.000 |
| perlin(L=4) | 0.40 | 433 | 0.540 | 0.140 | 0.581 | 0.258 | 1.000 | 0.000 |

## tetraglyme: 56 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 0.20 | 462 | 0.050 | 0.000 | 0.066 | 1.231 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.10 | 463 | 0.043 | 0.013 | 0.059 | 1.330 | 1.000 | 0.000 |
| perlin(L=8) | 0.20 | 484 | 0.077 | 0.016 | 0.102 | 1.108 | 1.000 | 0.000 |
| gauss-corr(L=8) | 0.10 | 488 | 0.077 | 0.016 | 0.102 | 1.030 | 1.000 | 0.000 |
| gauss-iid | 0.00 | 500 | 0.023 | 0.005 | 0.032 | 1.393 | 1.000 | 0.000 |

## ala3: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 1.00 | 94 | 0.975 | 0.550 | 0.988 | 0.013 | 1.000 | 0.000 |
| perlin(L=8) | 1.00 | 96 | 0.963 | 0.512 | 0.975 | 0.016 | 1.000 | 0.006 |
| gauss-iid | 1.00 | 98 | 0.963 | 0.487 | 0.988 | 0.014 | 1.000 | 0.000 |
| perlin(L=4) | 1.00 | 106 | 0.975 | 0.512 | 0.988 | 0.013 | 1.000 | 0.000 |
| gauss-corr(L=8) | 0.80 | 154 | 0.938 | 0.613 | 0.938 | 0.045 | 1.000 | 0.000 |
