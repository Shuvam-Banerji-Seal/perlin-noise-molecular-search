# Ablation and time-to-known-structure analysis

Targets are the exhaustive-grid minima within 3.0 kcal/mol of the grid's global minimum. Budget = 500 relaxations. Times are right-censored at the budget (counted as 501).

## Paired comparisons (pooled over molecule × seed)

Δrecall > 0 and Δt < 0 favour the *variant*. 95% CI: paired bootstrap. p: Wilcoxon; Holm within metric.

| question | variant − reference | Δ final recall [95% CI] | p (Holm) | Δ steps to known global min [95% CI] | p (Holm) |
|---|---|---|---|---|---|
| Perlin vs matched Gaussian (baseline) | perlin(L=4) − gauss-corr(L=4) | +0.030 [-0.018, +0.080] | 0.293 (1) | -2.0 [-47.8, +42.4] | 0.955 (1) |
| correlation vs none | gauss-corr(L=4) − gauss-iid | -0.033 [-0.084, +0.020] | 0.0224 (0.36) | +10.5 [-29.9, +51.6] | 0.905 (1) |
| Perlin vs i.i.d. | perlin(L=4) − gauss-iid | -0.004 [-0.054, +0.049] | 0.471 (1) | +8.6 [-32.7, +48.9] | 0.97 (1) |
| Perlin vs MCMM reset | perlin(L=4) − reset | +0.022 [-0.029, +0.080] | 0.337 (1) | -2.3 [-44.2, +38.9] | 0.855 (1) |
| marginals: Gaussianize Perlin | perlin-gaussmarg(L=4) − perlin(L=4) | +0.033 [+0.001, +0.067] | 0.152 (1) | -40.1 [-76.1, -4.6] | 0.0181 (0.31) |
| marginals: Perlin-shaped Gaussian | gauss-perlinmarg(L=4) − gauss-corr(L=4) | +0.036 [-0.016, +0.091] | 0.406 (1) | -33.5 [-72.3, +5.4] | 0.15 (1) |
| octaves: 1 vs 2 | perlin-oct1(L=4) − perlin(L=4) | +0.034 [-0.012, +0.083] | 0.624 (1) | +5.4 [-33.8, +45.7] | 0.855 (1) |
| octaves: 4 vs 2 | perlin-oct4(L=4) − perlin(L=4) | +0.003 [-0.042, +0.050] | 0.792 (1) | -14.1 [-50.2, +21.7] | 0.406 (1) |
| 1 octave vs its matched Gaussian | perlin-oct1(L=4) − gauss-corr-oct1(L=4) | +0.040 [-0.008, +0.088] | 0.333 (1) | +41.5 [-2.5, +85.2] | 0.0535 (0.86) |
| 4 octaves vs its matched Gaussian | perlin-oct4(L=4) − gauss-corr-oct4(L=4) | +0.009 [-0.037, +0.057] | 0.855 (1) | -3.4 [-46.0, +39.0] | 0.797 (1) |
| random rotation off | perlin-norot(L=4) − perlin(L=4) | +0.004 [-0.050, +0.059] | 0.846 (1) | -7.8 [-49.4, +34.0] | 0.722 (1) |
| no rotation vs its matched Gaussian | perlin-norot(L=4) − gauss-corr-norot(L=4) | +0.004 [-0.051, +0.058] | 0.796 (1) | +7.5 [-31.4, +46.3] | 0.932 (1) |
| topology: shuffled sites | perlin-shuffled(L=4) − perlin(L=4) | +0.011 [-0.032, +0.054] | 0.607 (1) | -30.2 [-71.2, +11.5] | 0.425 (1) |
| topology: shuffled (Gaussian) | gauss-corr-shuffled(L=4) − gauss-corr(L=4) | +0.039 [-0.013, +0.095] | 0.297 (1) | -14.2 [-52.5, +23.5] | 0.859 (1) |
| temporal coherence on | perlin-drift(L=4) − perlin(L=4) | -0.041 [-0.088, +0.008] | 0.00362 (0.062) | +3.3 [-39.4, +45.9] | 0.506 (1) |
| temporal coherence vs matched AR(1) Gaussian | perlin-drift(L=4) − gauss-ar1(L=4) | +0.005 [-0.041, +0.052] | 0.497 (1) | +1.2 [-40.0, +45.9] | 0.916 (1) |
| temporal coherence (Gaussian) | gauss-ar1(L=4) − gauss-corr(L=4) | -0.016 [-0.061, +0.029] | 0.227 (1) | +0.2 [-41.7, +40.4] | 0.728 (1) |

## ala1: 7 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin(L=4) | 1.00 | 4 | 0.914 | 0.721 | 1.000 | 0.006 | 1.000 | 0.000 |
| gauss-ar1(L=4) | 1.00 | 4 | 0.921 | 0.721 | 1.000 | 0.007 | 1.000 | 0.000 |
| gauss-perlinmarg(L=4) | 1.00 | 5 | 0.921 | 0.729 | 1.000 | 0.004 | 1.000 | 0.000 |
| perlin-oct4(L=4) | 1.00 | 6 | 0.936 | 0.757 | 1.000 | 0.005 | 1.000 | 0.000 |
| perlin-gaussmarg(L=4) | 1.00 | 6 | 0.957 | 0.771 | 1.000 | 0.002 | 1.000 | 0.000 |
| gauss-corr-shuffled(L=4) | 1.00 | 6 | 0.943 | 0.743 | 1.000 | 0.002 | 1.000 | 0.000 |
| gauss-corr-oct4(L=4) | 1.00 | 6 | 0.943 | 0.779 | 1.000 | 0.002 | 1.000 | 0.000 |
| gauss-corr-oct1(L=4) | 1.00 | 6 | 0.893 | 0.721 | 1.000 | 0.007 | 1.000 | 0.000 |
| gauss-corr(L=4) | 1.00 | 6 | 0.914 | 0.707 | 1.000 | 0.006 | 1.000 | 0.000 |
| gauss-iid | 1.00 | 7 | 0.993 | 0.814 | 1.000 | 0.002 | 1.000 | 0.000 |
| perlin-shuffled(L=4) | 1.00 | 7 | 0.914 | 0.714 | 1.000 | 0.005 | 1.000 | 0.000 |
| perlin-oct1(L=4) | 1.00 | 8 | 0.893 | 0.693 | 1.000 | 0.005 | 1.000 | 0.000 |
| gauss-corr-norot(L=4) | 1.00 | 8 | 0.929 | 0.764 | 1.000 | 0.004 | 1.000 | 0.000 |
| perlin-norot(L=4) | 1.00 | 9 | 0.900 | 0.714 | 1.000 | 0.007 | 1.000 | 0.000 |
| perlin-drift(L=4) | 1.00 | 10 | 0.921 | 0.729 | 1.000 | 0.006 | 1.000 | 0.000 |
| reset | 1.00 | 12 | 1.000 | 0.979 | 1.000 | 0.000 | 1.000 | 0.000 |

## decane: 193 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| gauss-corr-oct1(L=4) | 0.90 | 176 | 0.580 | 0.175 | 0.613 | 0.235 | 1.000 | 0.000 |
| perlin-norot(L=4) | 0.95 | 177 | 0.587 | 0.177 | 0.624 | 0.231 | 1.000 | 0.000 |
| perlin-gaussmarg(L=4) | 0.90 | 196 | 0.587 | 0.166 | 0.622 | 0.231 | 1.000 | 0.000 |
| perlin-oct4(L=4) | 0.80 | 207 | 0.602 | 0.181 | 0.639 | 0.222 | 1.000 | 0.000 |
| gauss-iid | 0.85 | 213 | 0.606 | 0.182 | 0.639 | 0.221 | 1.000 | 0.000 |
| gauss-perlinmarg(L=4) | 0.90 | 213 | 0.591 | 0.169 | 0.625 | 0.229 | 1.000 | 0.000 |
| perlin-shuffled(L=4) | 0.80 | 216 | 0.571 | 0.167 | 0.606 | 0.242 | 1.000 | 0.000 |
| gauss-corr-norot(L=4) | 0.90 | 226 | 0.580 | 0.165 | 0.615 | 0.236 | 1.000 | 0.000 |
| reset | 0.80 | 230 | 0.587 | 0.180 | 0.619 | 0.235 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.85 | 230 | 0.595 | 0.173 | 0.630 | 0.226 | 1.000 | 0.000 |
| perlin(L=4) | 0.90 | 256 | 0.602 | 0.183 | 0.633 | 0.223 | 1.000 | 0.000 |
| gauss-corr-oct4(L=4) | 0.65 | 258 | 0.582 | 0.171 | 0.618 | 0.235 | 1.000 | 0.000 |
| gauss-corr-shuffled(L=4) | 0.80 | 275 | 0.562 | 0.160 | 0.597 | 0.247 | 1.000 | 0.000 |
| gauss-ar1(L=4) | 0.70 | 324 | 0.503 | 0.135 | 0.543 | 0.282 | 1.000 | 0.000 |
| perlin-drift(L=4) | 0.55 | 350 | 0.480 | 0.129 | 0.523 | 0.295 | 1.000 | 0.000 |
| perlin-oct1(L=4) | 0.70 | 356 | 0.584 | 0.165 | 0.614 | 0.233 | 1.000 | 0.000 |

## tetraglyme: 56 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin-oct1(L=4) | 0.40 | 402 | 0.162 | 0.024 | 0.199 | 0.811 | 1.000 | 0.000 |
| gauss-ar1(L=4) | 0.25 | 404 | 0.100 | 0.026 | 0.135 | 1.033 | 1.000 | 0.000 |
| gauss-corr-oct4(L=4) | 0.40 | 418 | 0.140 | 0.023 | 0.172 | 0.906 | 1.000 | 0.000 |
| gauss-corr-oct1(L=4) | 0.30 | 419 | 0.163 | 0.033 | 0.189 | 0.794 | 1.000 | 0.000 |
| gauss-perlinmarg(L=4) | 0.45 | 423 | 0.122 | 0.023 | 0.156 | 0.936 | 1.000 | 0.000 |
| reset | 0.30 | 427 | 0.156 | 0.038 | 0.200 | 0.856 | 1.000 | 0.000 |
| gauss-corr-norot(L=4) | 0.35 | 430 | 0.163 | 0.032 | 0.210 | 0.768 | 1.000 | 0.000 |
| perlin(L=4) | 0.20 | 445 | 0.128 | 0.030 | 0.154 | 0.908 | 1.000 | 0.000 |
| perlin-oct4(L=4) | 0.20 | 445 | 0.096 | 0.021 | 0.120 | 1.032 | 1.000 | 0.000 |
| perlin-norot(L=4) | 0.30 | 453 | 0.151 | 0.025 | 0.194 | 0.827 | 1.000 | 0.000 |
| perlin-gaussmarg(L=4) | 0.20 | 455 | 0.114 | 0.021 | 0.139 | 0.904 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.30 | 464 | 0.117 | 0.011 | 0.145 | 0.913 | 1.000 | 0.000 |
| perlin-shuffled(L=4) | 0.15 | 483 | 0.096 | 0.015 | 0.132 | 0.926 | 1.000 | 0.000 |
| gauss-iid | 0.05 | 488 | 0.088 | 0.013 | 0.111 | 0.943 | 1.000 | 0.000 |
| gauss-corr-shuffled(L=4) | 0.10 | 488 | 0.109 | 0.017 | 0.133 | 0.976 | 1.000 | 0.000 |
| perlin-drift(L=4) | 0.05 | 498 | 0.052 | 0.011 | 0.063 | 1.088 | 1.000 | 0.000 |

## ala3: 8 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin-shuffled(L=4) | 0.95 | 104 | 0.906 | 0.525 | 0.938 | 0.121 | 1.000 | 0.000 |
| perlin-gaussmarg(L=4) | 0.95 | 142 | 0.900 | 0.556 | 0.963 | 0.067 | 1.000 | 0.000 |
| perlin-drift(L=4) | 0.90 | 153 | 0.819 | 0.419 | 0.887 | 0.187 | 1.000 | 0.000 |
| gauss-corr-oct1(L=4) | 0.95 | 157 | 0.856 | 0.537 | 0.906 | 0.152 | 1.000 | 0.000 |
| gauss-iid | 0.95 | 164 | 0.900 | 0.450 | 0.956 | 0.091 | 1.000 | 0.000 |
| gauss-corr-shuffled(L=4) | 0.85 | 168 | 0.856 | 0.438 | 0.906 | 0.168 | 1.000 | 0.000 |
| gauss-perlinmarg(L=4) | 0.85 | 180 | 0.794 | 0.406 | 0.856 | 0.254 | 1.000 | 0.000 |
| reset | 0.95 | 189 | 0.881 | 0.375 | 0.931 | 0.136 | 1.000 | 0.000 |
| gauss-corr-oct4(L=4) | 0.75 | 203 | 0.806 | 0.500 | 0.850 | 0.304 | 1.000 | 0.000 |
| gauss-corr-norot(L=4) | 0.75 | 217 | 0.800 | 0.438 | 0.838 | 0.317 | 1.000 | 0.000 |
| gauss-ar1(L=4) | 0.80 | 217 | 0.794 | 0.444 | 0.856 | 0.271 | 1.000 | 0.000 |
| gauss-corr(L=4) | 0.60 | 256 | 0.713 | 0.444 | 0.800 | 0.301 | 1.000 | 0.003 |
| perlin-oct4(L=4) | 0.85 | 259 | 0.844 | 0.431 | 0.931 | 0.115 | 1.000 | 0.000 |
| perlin(L=4) | 0.60 | 270 | 0.819 | 0.419 | 0.856 | 0.213 | 1.000 | 0.000 |
| perlin-norot(L=4) | 0.75 | 280 | 0.850 | 0.431 | 0.887 | 0.190 | 1.000 | 0.000 |
| perlin-oct1(L=4) | 0.80 | 282 | 0.825 | 0.450 | 0.869 | 0.240 | 1.000 | 0.000 |

## imatinib: 79 known low-energy minima

| method | P(find known global min) | restricted mean steps to it | final recall | recall @100 | COV-R@0.5Å | MAT-R (Å) | COV-P@0.5Å | MAT-P (Å) |
|---|---|---|---|---|---|---|---|---|
| perlin-oct1(L=4) | 0.90 | 217 | 0.719 | 0.097 | 0.797 | 0.315 | 0.941 | 0.037 |
| perlin-gaussmarg(L=4) | 0.75 | 239 | 0.623 | 0.118 | 0.694 | 0.590 | 0.945 | 0.034 |
| gauss-corr-shuffled(L=4) | 0.70 | 241 | 0.590 | 0.133 | 0.649 | 0.609 | 0.899 | 0.059 |
| perlin-drift(L=4) | 0.70 | 244 | 0.540 | 0.122 | 0.609 | 0.729 | 0.944 | 0.036 |
| perlin-oct4(L=4) | 0.65 | 250 | 0.554 | 0.115 | 0.609 | 0.778 | 0.912 | 0.052 |
| gauss-perlinmarg(L=4) | 0.75 | 259 | 0.620 | 0.118 | 0.681 | 0.595 | 0.950 | 0.031 |
| perlin(L=4) | 0.65 | 263 | 0.553 | 0.146 | 0.614 | 0.777 | 0.905 | 0.056 |
| perlin-shuffled(L=4) | 0.70 | 277 | 0.582 | 0.075 | 0.648 | 0.690 | 0.947 | 0.033 |
| perlin-norot(L=4) | 0.65 | 280 | 0.547 | 0.087 | 0.613 | 0.780 | 0.906 | 0.055 |
| gauss-corr-norot(L=4) | 0.65 | 281 | 0.544 | 0.054 | 0.604 | 0.801 | 0.952 | 0.030 |
| gauss-corr(L=4) | 0.65 | 291 | 0.528 | 0.104 | 0.592 | 0.825 | 0.966 | 0.022 |
| gauss-corr-oct4(L=4) | 0.65 | 299 | 0.513 | 0.098 | 0.572 | 0.812 | 0.954 | 0.029 |
| gauss-corr-oct1(L=4) | 0.60 | 300 | 0.491 | 0.081 | 0.556 | 0.894 | 0.958 | 0.026 |
| gauss-ar1(L=4) | 0.60 | 300 | 0.466 | 0.090 | 0.530 | 0.921 | 0.954 | 0.029 |
| gauss-iid | 0.55 | 325 | 0.445 | 0.048 | 0.503 | 0.989 | 0.859 | 0.081 |
| reset | 0.40 | 392 | 0.281 | 0.038 | 0.323 | 0.987 | 0.750 | 0.140 |
