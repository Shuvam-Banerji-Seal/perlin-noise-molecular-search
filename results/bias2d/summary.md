# Perlin bias on the alanine-dipeptide (phi, psi) surface (formulation A1)

exact dF(phi>0 - phi<0) = 2.402 kcal/mol; min E on phi=0: 13.43, on phi=180: 2.86 kcal/mol; kT = 0.596

| condition | transitions / 1e5 steps | walkers crossed | dF mean | dF RMSE | finite dF |
|---|---|---|---|---|---|
| none@0 | 0.02 ± 0.02 | 0.06 | 3.84 | 1.59 | 1.00 |
| metad@0 | 64.89 ± 1.12 | 1.00 | 2.42 | 0.12 | 1.00 |
| perlin@2 | 0.00 ± 0.00 | 0.00 | 4.12 | 1.73 | 0.75 |
| perlin@4 | 0.00 ± 0.00 | 0.00 | 5.85 | 3.77 | 0.69 |
| perlin@6 | 0.00 ± 0.00 | 0.00 | 8.18 | 6.48 | 0.47 |
| grf@2 | 0.05 ± 0.03 | 0.09 | 4.14 | 2.07 | 0.88 |
| grf@4 | 0.00 ± 0.00 | 0.00 | 9.16 | 7.17 | 0.50 |
| grf@6 | 0.02 ± 0.01 | 0.06 | 11.70 | 10.15 | 0.34 |
| perlin-drift@2 | 0.02 ± 0.02 | 0.06 | 4.06 | 1.95 | 1.00 |
| perlin-drift@4 | 0.41 ± 0.13 | 0.31 | 3.24 | 2.66 | 1.00 |
| perlin-drift@6 | 8.01 ± 0.97 | 0.78 | 1.02 | 3.93 | 1.00 |
| grf-drift@2 | 0.29 ± 0.07 | 0.44 | 2.98 | 1.32 | 1.00 |
| grf-drift@4 | 3.66 ± 0.20 | 1.00 | 1.75 | 1.50 | 1.00 |
| grf-drift@6 | 8.02 ± 0.31 | 1.00 | 2.07 | 1.61 | 1.00 |
