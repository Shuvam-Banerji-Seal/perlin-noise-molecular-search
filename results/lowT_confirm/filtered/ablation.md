# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 500 relaxations. Times are right-censored at the budget (counted as 501).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| low-T confirmation | perlin(L=4) − gauss-corr(L=4) | -0.008 [-0.039, +0.025] | 0.453 (0.45) | +27.9 [-21.8, +78.6] | 0.364 (0.36) |

## decane: 193 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 1.00 | 105 | 0.608 | 0.209 | 0.628 | 0.220 | 1.000 | 0.000 |
| perlin(L=4) | 0.97 | 109 | 0.615 | 0.223 | 0.638 | 0.217 | 1.000 | 0.000 |

## tetraglyme: 56 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.70 | 323 | 0.217 | 0.042 | 0.253 | 0.867 | 1.000 | 0.000 |
| perlin(L=4) | 0.37 | 402 | 0.162 | 0.023 | 0.192 | 0.906 | 1.000 | 0.000 |

## ala3: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 0.60 | 258 | 0.658 | 0.471 | 0.721 | 0.512 | 1.000 | 0.002 |
| perlin(L=4) | 0.63 | 258 | 0.683 | 0.421 | 0.767 | 0.446 | 1.000 | 0.000 |
