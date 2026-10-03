# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 400 relaxations. Times are right-censored at the budget (counted as 401).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Cartesian Perlin vs matched field (L=2 A) | cart-perlin@2 − cart-gauss-field@2 | -0.005 [-0.047, +0.037] | 0.707 (1) | +9.8 [-34.5, +57.0] | 0.765 (1) |
| Cartesian Perlin vs matched field (L=4 A) | cart-perlin@4 − cart-gauss-field@4 | +0.009 [-0.033, +0.053] | 0.866 (1) | -10.2 [-64.5, +45.3] | 0.79 (1) |
| Cartesian Perlin (L=2 A) vs i.i.d. | cart-perlin@2 − cart-iid | -0.002 [-0.056, +0.060] | 0.396 (1) | +10.8 [-39.1, +60.0] | 0.627 (1) |
| Cartesian correlated field vs i.i.d. | cart-gauss-field@2 − cart-iid | +0.003 [-0.049, +0.059] | 0.793 (1) | +0.9 [-58.0, +60.8] | 0.844 (1) |

## imatinib: 79 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| cart-iid | 0.90 | 164 | 0.377 | 0.115 | 0.447 | 0.537 | 0.954 | 0.027 |
| cart-perlin@2 | 0.80 | 206 | 0.375 | 0.127 | 0.432 | 0.575 | 0.916 | 0.050 |
| cart-gauss-field@2 | 0.80 | 247 | 0.389 | 0.099 | 0.444 | 0.536 | 0.936 | 0.039 |
| cart-gauss-field@4 | 0.70 | 272 | 0.319 | 0.082 | 0.370 | 0.774 | 0.946 | 0.032 |
| cart-perlin@4 | 0.50 | 279 | 0.329 | 0.082 | 0.392 | 0.707 | 0.938 | 0.037 |

## ala3: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| cart-gauss-field@2 | 0.60 | 195 | 0.550 | 0.312 | 0.550 | 0.876 | 1.000 | 0.044 |
| cart-perlin@2 | 0.60 | 229 | 0.550 | 0.263 | 0.575 | 0.872 | 1.000 | 0.042 |
| cart-iid | 0.50 | 262 | 0.537 | 0.325 | 0.550 | 0.955 | 1.000 | 0.033 |
| cart-perlin@4 | 0.40 | 308 | 0.412 | 0.250 | 0.425 | 1.262 | 1.000 | 0.031 |
| cart-gauss-field@4 | 0.20 | 350 | 0.350 | 0.175 | 0.350 | 1.509 | 1.000 | 0.009 |

## ala1: 7 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| cart-iid | 1.00 | 4 | 0.843 | 0.657 | 0.914 | 0.095 | 1.000 | 0.000 |
| cart-perlin@2 | 1.00 | 6 | 0.829 | 0.643 | 0.914 | 0.117 | 1.000 | 0.000 |
| cart-gauss-field@2 | 1.00 | 11 | 0.829 | 0.614 | 0.943 | 0.088 | 1.000 | 0.001 |
| cart-gauss-field@4 | 1.00 | 11 | 0.729 | 0.500 | 0.800 | 0.242 | 1.000 | 0.000 |
| cart-perlin@4 | 1.00 | 20 | 0.714 | 0.571 | 0.771 | 0.259 | 1.000 | 0.000 |

## decane: 193 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| cart-gauss-field@2 | 0.60 | 230 | 0.256 | 0.089 | 0.290 | 0.466 | 1.000 | 0.000 |
| cart-iid | 0.60 | 249 | 0.254 | 0.080 | 0.290 | 0.464 | 1.000 | 0.000 |
| cart-perlin@2 | 0.60 | 281 | 0.252 | 0.073 | 0.288 | 0.466 | 1.000 | 0.000 |
| cart-perlin@4 | 0.20 | 349 | 0.115 | 0.028 | 0.132 | 0.607 | 1.000 | 0.000 |
| cart-gauss-field@4 | 0.10 | 362 | 0.138 | 0.046 | 0.164 | 0.581 | 1.000 | 0.000 |
