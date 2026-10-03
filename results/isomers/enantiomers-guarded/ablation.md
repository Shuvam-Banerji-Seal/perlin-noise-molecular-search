# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 300 relaxations. Times are right-censored at the budget (counted as 301).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched Gaussian (baseline) | perlin(L=4) − gauss-corr(L=4) | -0.028 [-0.069, +0.014] | 0.158 (0.24) | +29.6 [+2.0, +57.3] | 0.312 (0.62) |
| correlation vs none | gauss-corr(L=4) − gauss-iid | -0.031 [-0.075, +0.012] | 0.12 (0.24) | +3.8 [-16.5, +24.5] | 0.434 (0.62) |
| Perlin vs i.i.d. | perlin(L=4) − gauss-iid | -0.059 [-0.104, -0.014] | 0.00988 (0.03) | +33.4 [+9.8, +58.1] | 0.0363 (0.15) |
| Perlin vs MCMM reset | perlin(L=4) − reset | -0.065 [-0.105, -0.025] | 0.000938 (0.0038) | +21.2 [-6.0, +50.4] | 0.193 (0.58) |

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
| perlin(L=4) | 0.80 | 129 | 0.552 | 0.370 | 0.642 | 0.301 | 1.000 | 0.005 |
| reset | 0.90 | 138 | 0.568 | 0.347 | 0.648 | 0.286 | 1.000 | 0.005 |
| gauss-corr(L=4) | 0.80 | 141 | 0.532 | 0.300 | 0.660 | 0.282 | 1.000 | 0.003 |
| gauss-iid | 0.60 | 165 | 0.645 | 0.395 | 0.750 | 0.200 | 1.000 | 0.002 |

## propranolol-S: 40 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-iid | 1.00 | 43 | 0.565 | 0.318 | 0.652 | 0.274 | 1.000 | 0.005 |
| gauss-corr(L=4) | 0.90 | 94 | 0.573 | 0.320 | 0.663 | 0.276 | 1.000 | 0.005 |
| reset | 0.70 | 143 | 0.595 | 0.295 | 0.690 | 0.255 | 1.000 | 0.006 |
| perlin(L=4) | 0.60 | 149 | 0.467 | 0.240 | 0.570 | 0.368 | 1.000 | 0.006 |

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
