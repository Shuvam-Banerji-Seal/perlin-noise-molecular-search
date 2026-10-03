# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 500 relaxations. Times are right-censored at the budget (counted as 501).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched (L=4), kT=0.3 | perlin(L=4) − gauss-corr(L=4) | +0.050 [+0.014, +0.091] | 0.0263 (0.072) | -6.0 [-88.8, +78.0] | 0.949 (1) |
| Perlin vs matched (L=8), kT=0.3 | perlin(L=8) − gauss-corr(L=8) | +0.007 [-0.041, +0.058] | 0.872 (0.87) | -8.4 [-90.7, +70.6] | 0.889 (1) |
| correlation vs none, kT=0.3 | gauss-corr(L=8) − gauss-iid | -0.048 [-0.086, -0.011] | 0.0238 (0.072) | -51.9 [-119.3, +18.5] | 0.162 (0.48) |

## decane: 193 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=4) | 1.00 | 70 | 0.614 | 0.225 | 0.635 | 0.217 | 1.000 | 0.000 |
| gauss-iid | 1.00 | 71 | 0.615 | 0.217 | 0.631 | 0.221 | 1.000 | 0.000 |
| perlin(L=8) | 1.00 | 95 | 0.460 | 0.177 | 0.482 | 0.313 | 1.000 | 0.000 |
| perlin(L=4) | 1.00 | 133 | 0.647 | 0.218 | 0.668 | 0.197 | 1.000 | 0.000 |
| gauss-corr(L=8) | 1.00 | 139 | 0.477 | 0.172 | 0.499 | 0.302 | 1.000 | 0.000 |

## tetraglyme: 56 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr(L=8) | 1.00 | 228 | 0.227 | 0.037 | 0.259 | 0.896 | 1.000 | 0.000 |
| perlin(L=4) | 0.80 | 300 | 0.229 | 0.020 | 0.266 | 0.885 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.60 | 333 | 0.211 | 0.032 | 0.237 | 0.892 | 1.000 | 0.000 |
| perlin(L=8) | 0.70 | 344 | 0.189 | 0.032 | 0.225 | 0.821 | 1.000 | 0.000 |
| gauss-iid | 0.30 | 461 | 0.184 | 0.016 | 0.227 | 0.775 | 1.000 | 0.000 |

## ala3: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=8) | 0.50 | 266 | 0.688 | 0.450 | 0.800 | 0.316 | 1.000 | 0.000 |
| perlin(L=4) | 0.50 | 269 | 0.650 | 0.425 | 0.750 | 0.435 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.40 | 318 | 0.550 | 0.388 | 0.675 | 0.558 | 1.000 | 0.000 |
| gauss-iid | 0.40 | 353 | 0.662 | 0.375 | 0.775 | 0.338 | 1.000 | 0.000 |
| gauss-corr(L=8) | 0.30 | 363 | 0.613 | 0.450 | 0.675 | 0.564 | 1.000 | 0.000 |
