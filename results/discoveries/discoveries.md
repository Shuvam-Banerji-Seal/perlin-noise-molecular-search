# Discoveries: pooled minima across all runs and the grid reference

## decane

- unique minima (all sources): 1808; within 3.0 kcal/mol of E* = -6.299: **193**
- grid reference low-energy minima: 193; missed by every method: **0**
- low-energy minima found by basin hopping but **not** by the grid: **0** (ΔE: [])
- global minimum in grid reference: True
- global minimum found by: etkdg @60, gauss-corr(L=2) @60, gauss-corr(L=4) @120, gauss-corr(L=4) @30, gauss-corr(L=4) @60, gauss-corr(L=8) @120, gauss-corr(L=8) @30, gauss-corr(L=8) @60, gauss-iid @120, gauss-iid @30, gauss-iid @60, grid, perlin(L=2) @60, perlin(L=4) @120, perlin(L=4) @30, perlin(L=4) @60, perlin(L=8) @120, perlin(L=8) @30, perlin(L=8) @60, reset @60
- H-rebuild energy check: max |ΔE| = 0.000 kcal/mol

| method @ sigma | coverage (common ref) | recall vs grid |
|---|---|---|
| etkdg @60 | 0.716 ± 0.022 | 0.716 ± 0.022 |
| gauss-corr(L=2) @60 | 0.616 ± 0.025 | 0.616 ± 0.025 |
| perlin(L=4) @60 | 0.616 ± 0.021 | 0.616 ± 0.021 |
| perlin(L=2) @60 | 0.612 ± 0.028 | 0.612 ± 0.028 |
| gauss-iid @60 | 0.607 ± 0.041 | 0.607 ± 0.041 |
| gauss-corr(L=4) @60 | 0.597 ± 0.031 | 0.597 ± 0.031 |
| reset @60 | 0.575 ± 0.047 | 0.575 ± 0.047 |
| perlin(L=8) @60 | 0.555 ± 0.055 | 0.555 ± 0.055 |
| gauss-corr(L=8) @60 | 0.541 ± 0.032 | 0.541 ± 0.032 |
| gauss-corr(L=8) @120 | 0.507 ± 0.026 | 0.507 ± 0.026 |
| perlin(L=4) @120 | 0.493 ± 0.037 | 0.493 ± 0.037 |
| gauss-corr(L=4) @120 | 0.483 ± 0.046 | 0.483 ± 0.046 |
| perlin(L=8) @120 | 0.482 ± 0.028 | 0.482 ± 0.028 |
| gauss-iid @120 | 0.456 ± 0.031 | 0.456 ± 0.031 |
| gauss-iid @30 | 0.351 ± 0.045 | 0.351 ± 0.045 |
| perlin(L=4) @30 | 0.319 ± 0.031 | 0.319 ± 0.031 |
| gauss-corr(L=4) @30 | 0.315 ± 0.022 | 0.315 ± 0.022 |
| perlin(L=8) @30 | 0.313 ± 0.034 | 0.313 ± 0.034 |
| gauss-corr(L=8) @30 | 0.294 ± 0.044 | 0.294 ± 0.044 |

## tetraglyme

- unique minima (all sources): 164924; within 3.0 kcal/mol of E* = 88.954: **56**
- grid reference low-energy minima: 56; missed by every method: **0**
- low-energy minima found by basin hopping but **not** by the grid: **0** (ΔE: [])
- global minimum in grid reference: True
- global minimum found by: gauss-corr(L=2) @60, gauss-corr(L=4) @30, gauss-corr(L=4) @60, gauss-corr(L=8) @30, gauss-corr(L=8) @60, gauss-iid @30, gauss-iid @60, grid, perlin(L=4) @30, perlin(L=4) @60, perlin(L=8) @30, perlin(L=8) @60, reset @60
- H-rebuild energy check: max |ΔE| = 0.000 kcal/mol

| method @ sigma | coverage (common ref) | recall vs grid |
|---|---|---|
| gauss-iid @30 | 0.421 ± 0.039 | 0.421 ± 0.039 |
| gauss-corr(L=4) @30 | 0.325 ± 0.036 | 0.325 ± 0.036 |
| perlin(L=4) @30 | 0.289 ± 0.063 | 0.289 ± 0.063 |
| perlin(L=8) @30 | 0.266 ± 0.051 | 0.266 ± 0.051 |
| gauss-corr(L=8) @30 | 0.257 ± 0.053 | 0.257 ± 0.053 |
| gauss-corr(L=8) @60 | 0.205 ± 0.065 | 0.205 ± 0.065 |
| perlin(L=8) @60 | 0.188 ± 0.044 | 0.188 ± 0.044 |
| reset @60 | 0.155 ± 0.116 | 0.155 ± 0.116 |
| perlin(L=4) @60 | 0.146 ± 0.059 | 0.146 ± 0.059 |
| gauss-corr(L=4) @60 | 0.120 ± 0.100 | 0.120 ± 0.100 |
| gauss-iid @60 | 0.089 ± 0.039 | 0.089 ± 0.039 |
| gauss-corr(L=2) @60 | 0.089 ± 0.061 | 0.089 ± 0.061 |
| perlin(L=2) @60 | 0.089 ± 0.024 | 0.089 ± 0.024 |
| gauss-corr(L=8) @120 | 0.025 ± 0.026 | 0.025 ± 0.026 |
| perlin(L=8) @120 | 0.020 ± 0.018 | 0.020 ± 0.018 |
| gauss-corr(L=4) @120 | 0.005 ± 0.012 | 0.005 ± 0.012 |
| gauss-iid @120 | 0.004 ± 0.008 | 0.004 ± 0.008 |
| perlin(L=4) @120 | 0.004 ± 0.008 | 0.004 ± 0.008 |
| etkdg @60 | 0.000 ± 0.000 | 0.000 ± 0.000 |

## ala3

- unique minima (all sources): 992; within 3.0 kcal/mol of E* = -6.459: **9**
- grid reference low-energy minima: 9; missed by every method: **1**
- low-energy minima found by basin hopping but **not** by the grid: **0** (ΔE: [])
- global minimum in grid reference: True
- global minimum found by: gauss-corr(L=2) @60, gauss-corr(L=4) @120, gauss-corr(L=4) @30, gauss-corr(L=4) @60, gauss-corr(L=8) @120, gauss-corr(L=8) @30, gauss-corr(L=8) @60, gauss-iid @120, gauss-iid @30, gauss-iid @60, grid, perlin(L=2) @60, perlin(L=4) @120, perlin(L=4) @30, perlin(L=4) @60, perlin(L=8) @120, perlin(L=8) @30, perlin(L=8) @60, reset @60
- H-rebuild energy check: max |ΔE| = 0.000 kcal/mol

| method @ sigma | coverage (common ref) | recall vs grid |
|---|---|---|
| gauss-corr(L=2) @60 | 0.856 ± 0.075 | 0.856 ± 0.075 |
| gauss-iid @60 | 0.822 ± 0.094 | 0.822 ± 0.094 |
| gauss-corr(L=4) @120 | 0.800 ± 0.088 | 0.800 ± 0.088 |
| gauss-corr(L=8) @120 | 0.800 ± 0.070 | 0.800 ± 0.070 |
| gauss-iid @120 | 0.789 ± 0.097 | 0.789 ± 0.097 |
| gauss-corr(L=8) @60 | 0.778 ± 0.128 | 0.778 ± 0.128 |
| perlin(L=4) @120 | 0.767 ± 0.082 | 0.767 ± 0.082 |
| reset @60 | 0.744 ± 0.158 | 0.744 ± 0.158 |
| perlin(L=2) @60 | 0.744 ± 0.158 | 0.744 ± 0.158 |
| perlin(L=8) @60 | 0.744 ± 0.158 | 0.744 ± 0.158 |
| perlin(L=8) @120 | 0.700 ± 0.091 | 0.700 ± 0.091 |
| perlin(L=4) @60 | 0.678 ± 0.212 | 0.678 ± 0.212 |
| gauss-corr(L=4) @60 | 0.633 ± 0.189 | 0.633 ± 0.189 |
| gauss-iid @30 | 0.433 ± 0.169 | 0.433 ± 0.169 |
| gauss-corr(L=8) @30 | 0.400 ± 0.223 | 0.400 ± 0.223 |
| perlin(L=8) @30 | 0.356 ± 0.187 | 0.356 ± 0.187 |
| gauss-corr(L=4) @30 | 0.344 ± 0.169 | 0.344 ± 0.169 |
| perlin(L=4) @30 | 0.333 ± 0.157 | 0.333 ± 0.157 |
| etkdg @60 | 0.000 ± 0.000 | 0.000 ± 0.000 |

## imatinib

- unique minima (all sources): 5030; within 3.0 kcal/mol of E* = 64.812: **88**
- grid reference low-energy minima: 79; missed by every method: **0**
- low-energy minima found by basin hopping but **not** by the grid: **9** (ΔE: [1.575, 2.56, 2.56, 2.609, 2.769, 2.769, 2.879, 2.882, 2.952])
- global minimum in grid reference: True
- global minimum found by: etkdg @60, gauss-corr(L=2) @60, gauss-corr(L=4) @120, gauss-corr(L=4) @30, gauss-corr(L=4) @60, gauss-corr(L=8) @120, gauss-corr(L=8) @30, gauss-corr(L=8) @60, gauss-iid @120, gauss-iid @60, grid, perlin(L=2) @60, perlin(L=4) @120, perlin(L=4) @30, perlin(L=4) @60, perlin(L=8) @120, perlin(L=8) @30, perlin(L=8) @60, reset @60
- H-rebuild energy check: max |ΔE| = 0.000 kcal/mol

| method @ sigma | coverage (common ref) | recall vs grid |
|---|---|---|
| gauss-corr(L=8) @60 | 0.582 ± 0.389 | 0.609 ± 0.406 |
| gauss-corr(L=2) @60 | 0.573 ± 0.383 | 0.603 ± 0.401 |
| gauss-corr(L=4) @120 | 0.533 ± 0.356 | 0.553 ± 0.369 |
| perlin(L=4) @120 | 0.532 ± 0.354 | 0.549 ± 0.365 |
| gauss-iid @60 | 0.481 ± 0.395 | 0.491 ± 0.408 |
| perlin(L=8) @60 | 0.476 ± 0.393 | 0.491 ± 0.409 |
| perlin(L=2) @60 | 0.444 ± 0.372 | 0.463 ± 0.386 |
| perlin(L=4) @60 | 0.440 ± 0.440 | 0.452 ± 0.456 |
| perlin(L=8) @120 | 0.434 ± 0.363 | 0.457 ± 0.381 |
| gauss-corr(L=8) @120 | 0.387 ± 0.386 | 0.404 ± 0.401 |
| gauss-corr(L=4) @60 | 0.349 ± 0.351 | 0.371 ± 0.372 |
| gauss-iid @120 | 0.306 ± 0.367 | 0.314 ± 0.374 |
| reset @60 | 0.268 ± 0.365 | 0.275 ± 0.377 |
| gauss-corr(L=8) @30 | 0.233 ± 0.266 | 0.242 ± 0.283 |
| perlin(L=4) @30 | 0.222 ± 0.281 | 0.234 ± 0.308 |
| perlin(L=8) @30 | 0.152 ± 0.186 | 0.162 ± 0.201 |
| gauss-corr(L=4) @30 | 0.150 ± 0.231 | 0.159 ± 0.248 |
| etkdg @60 | 0.112 ± 0.022 | 0.124 ± 0.022 |
| gauss-iid @30 | 0.033 ± 0.082 | 0.037 ± 0.091 |
