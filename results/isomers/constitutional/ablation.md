# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 300 relaxations. Times are right-censored at the budget (counted as 301).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched Gaussian (baseline) | perlin(L=4) − gauss-corr(L=4) | +0.001 [-0.024, +0.023] | 0.375 (1) | +5.7 [-18.0, +31.0] | 0.448 (1) |
| correlation vs none | gauss-corr(L=4) − gauss-iid | -0.017 [-0.034, -0.001] | 0.0555 (0.22) | -7.2 [-31.3, +17.8] | 0.605 (1) |
| Perlin vs i.i.d. | perlin(L=4) − gauss-iid | -0.015 [-0.037, +0.003] | 0.494 (1) | -1.4 [-30.1, +28.6] | 0.489 (1) |
| Perlin vs MCMM reset | perlin(L=4) − reset | +0.012 [-0.016, +0.040] | 0.464 (1) | +12.7 [-9.8, +36.7] | 0.44 (1) |

## n-decane: 193 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-iid | 0.80 | 159 | 0.443 | 0.180 | 0.479 | 0.320 | 1.000 | 0.000 |
| reset | 0.70 | 187 | 0.427 | 0.184 | 0.467 | 0.335 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.60 | 194 | 0.433 | 0.174 | 0.475 | 0.325 | 1.000 | 0.000 |
| perlin(L=4) | 0.70 | 231 | 0.440 | 0.191 | 0.481 | 0.320 | 1.000 | 0.000 |

## 2-methylnonane: 179 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 1.00 | 90 | 0.512 | 0.230 | 0.542 | 0.281 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 94 | 0.494 | 0.207 | 0.517 | 0.292 | 1.000 | 0.000 |
| reset | 1.00 | 94 | 0.498 | 0.204 | 0.521 | 0.296 | 1.000 | 0.000 |
| gauss-iid | 0.80 | 108 | 0.506 | 0.216 | 0.527 | 0.287 | 1.000 | 0.000 |

## 2,7-dimethyloctane: 39 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| reset | 1.00 | 22 | 0.928 | 0.662 | 0.936 | 0.044 | 1.000 | 0.000 |
| perlin(L=4) | 1.00 | 41 | 0.905 | 0.592 | 0.921 | 0.054 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 49 | 0.918 | 0.585 | 0.938 | 0.046 | 1.000 | 0.000 |
| gauss-iid | 1.00 | 57 | 0.928 | 0.638 | 0.941 | 0.041 | 1.000 | 0.000 |

## 4-propylheptane: 25 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 1.00 | 42 | 0.876 | 0.532 | 0.896 | 0.074 | 1.000 | 0.005 |
| gauss-corr(L=4) | 1.00 | 42 | 0.848 | 0.552 | 0.872 | 0.088 | 1.000 | 0.005 |
| reset | 1.00 | 59 | 0.852 | 0.540 | 0.904 | 0.082 | 1.000 | 0.009 |
| gauss-iid | 1.00 | 80 | 0.900 | 0.572 | 0.932 | 0.057 | 1.000 | 0.005 |

## 3,3-diethylhexane: 3 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| reset | 1.00 | 21 | 0.933 | 0.767 | 0.933 | 0.058 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 39 | 1.000 | 0.933 | 1.000 | 0.000 | 1.000 | 0.000 |
| perlin(L=4) | 1.00 | 43 | 0.967 | 0.833 | 0.967 | 0.029 | 1.000 | 0.000 |
| gauss-iid | 1.00 | 50 | 1.000 | 0.833 | 1.000 | 0.000 | 1.000 | 0.000 |
