# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 300 relaxations. Times are right-censored at the budget (counted as 301).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched Gaussian (baseline) | perlin(L=4) − gauss-corr(L=4) | -0.040 [-0.084, +0.004] | 0.0954 (0.19) | +28.8 [+3.0, +55.2] | 0.2 (0.6) |
| correlation vs none | gauss-corr(L=4) − gauss-iid | -0.012 [-0.060, +0.035] | 0.409 (0.41) | +2.0 [-18.6, +21.5] | 0.499 (0.6) |
| Perlin vs i.i.d. | perlin(L=4) − gauss-iid | -0.052 [-0.099, -0.006] | 0.0319 (0.096) | +30.8 [+9.1, +54.0] | 0.0124 (0.049) |
| Perlin vs MCMM reset | perlin(L=4) − reset | -0.067 [-0.108, -0.028] | 0.000576 (0.0023) | +19.9 [-7.3, +48.7] | 0.254 (0.6) |

## ala1-S: 7 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 1.00 | 4 | 0.886 | 0.729 | 1.000 | 0.007 | 1.000 | 0.000 |
| gauss-iid | 1.00 | 5 | 0.971 | 0.800 | 1.000 | 0.007 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 6 | 0.886 | 0.714 | 1.000 | 0.010 | 1.000 | 0.000 |
| reset | 1.00 | 15 | 1.000 | 0.986 | 1.000 | 0.000 | 1.000 | 0.000 |

## ala1-R: 6 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-iid | 1.00 | 5 | 0.983 | 0.833 | 1.000 | 0.004 | 1.000 | 0.005 |
| gauss-corr(L=4) | 1.00 | 8 | 0.917 | 0.800 | 0.983 | 0.026 | 1.000 | 0.001 |
| perlin(L=4) | 1.00 | 8 | 1.000 | 0.783 | 1.000 | 0.000 | 1.000 | 0.002 |
| reset | 1.00 | 16 | 1.000 | 0.983 | 1.000 | 0.000 | 1.000 | 0.005 |

## ibuprofen-R: 18 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-iid | 1.00 | 5 | 0.994 | 0.922 | 0.994 | 0.003 | 1.000 | 0.000 |
| reset | 1.00 | 9 | 1.000 | 0.994 | 1.000 | 0.000 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 10 | 1.000 | 0.928 | 1.000 | 0.000 | 1.000 | 0.000 |
| perlin(L=4) | 1.00 | 12 | 1.000 | 0.939 | 1.000 | 0.000 | 1.000 | 0.000 |

## ibuprofen-S: 22 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| reset | 1.00 | 3 | 0.818 | 0.805 | 1.000 | 0.063 | 1.000 | 0.000 |
| perlin(L=4) | 1.00 | 4 | 0.818 | 0.773 | 1.000 | 0.063 | 1.000 | 0.000 |
| gauss-iid | 1.00 | 8 | 0.814 | 0.755 | 0.995 | 0.065 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 14 | 0.818 | 0.759 | 1.000 | 0.063 | 1.000 | 0.000 |

## propranolol-R: 40 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.80 | 133 | 0.540 | 0.300 | 0.655 | 0.283 | 0.978 | 0.024 |
| reset | 0.80 | 134 | 0.518 | 0.352 | 0.603 | 0.317 | 0.974 | 0.026 |
| perlin(L=4) | 0.70 | 139 | 0.505 | 0.332 | 0.610 | 0.332 | 0.926 | 0.067 |
| gauss-iid | 0.70 | 149 | 0.582 | 0.370 | 0.682 | 0.245 | 0.938 | 0.051 |

## propranolol-S: 40 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-iid | 0.90 | 64 | 0.500 | 0.263 | 0.600 | 0.318 | 0.921 | 0.069 |
| gauss-corr(L=4) | 0.90 | 94 | 0.588 | 0.320 | 0.677 | 0.266 | 1.000 | 0.005 |
| perlin(L=4) | 0.80 | 124 | 0.445 | 0.250 | 0.555 | 0.372 | 0.993 | 0.012 |
| reset | 0.70 | 143 | 0.595 | 0.295 | 0.690 | 0.254 | 0.968 | 0.031 |

## ala3-SSS: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-iid | 0.80 | 141 | 0.775 | 0.438 | 0.838 | 0.301 | 1.000 | 0.000 |
| reset | 0.80 | 156 | 0.738 | 0.388 | 0.812 | 0.346 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.60 | 174 | 0.662 | 0.438 | 0.762 | 0.412 | 1.000 | 0.000 |
| perlin(L=4) | 0.20 | 251 | 0.600 | 0.362 | 0.650 | 0.608 | 1.000 | 0.000 |

## ala3-RRR: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.90 | 76 | 0.750 | 0.512 | 0.850 | 0.325 | 1.000 | 0.000 |
| reset | 0.90 | 110 | 0.713 | 0.375 | 0.787 | 0.456 | 1.000 | 0.000 |
| gauss-iid | 0.80 | 121 | 0.637 | 0.425 | 0.725 | 0.538 | 1.000 | 0.000 |
| perlin(L=4) | 0.40 | 203 | 0.588 | 0.463 | 0.650 | 0.637 | 1.000 | 0.000 |
