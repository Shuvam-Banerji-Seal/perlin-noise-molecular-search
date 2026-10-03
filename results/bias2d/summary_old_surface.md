# Perlin bias on the alanine-dipeptide (phi, psi) surface (formulation A1)

exact dF(phi>0 - phi<0) = 1.064 kcal/mol; min E on phi=0: 13.43, on phi=180: 1.43 kcal/mol; kT = 0.596

| condition | transitions / 1e5 steps | walkers crossed | dF mean | dF RMSE | finite dF |
|---|---|---|---|---|---|
| none@0 | 0.00 ± 0.00 | 0.00 | 4.31 | 3.26 | 1.00 |
| metad@0 | 18.48 ± 0.47 | 1.00 | 0.94 | 0.25 | 1.00 |
| perlin@2 | 0.00 ± 0.00 | 0.00 | 4.35 | 3.43 | 0.66 |
| perlin@4 | 0.00 ± 0.00 | 0.00 | 7.18 | 6.31 | 0.25 |
| perlin@6 | 0.00 ± 0.00 | 0.00 | 12.73 | 12.01 | 0.19 |
| grf@2 | 0.00 ± 0.00 | 0.00 | 5.40 | 4.67 | 0.34 |
| grf@4 | 0.00 ± 0.00 | 0.00 | 11.91 | 10.86 | 0.25 |
| grf@6 | 0.00 ± 0.00 | 0.00 | 16.98 | 16.00 | 0.19 |
| perlin-drift@2 | 0.06 ± 0.03 | 0.16 | 3.84 | 3.24 | 1.00 |
| perlin-drift@4 | 1.70 ± 0.37 | 0.66 | 1.29 | 3.38 | 1.00 |
| perlin-drift@6 | 6.59 ± 1.09 | 0.94 | 1.22 | 4.03 | 1.00 |
| grf-drift@2 | 0.29 ± 0.07 | 0.41 | 2.28 | 2.46 | 1.00 |
| grf-drift@4 | 2.72 ± 0.17 | 1.00 | 1.53 | 1.59 | 1.00 |
| grf-drift@6 | 6.61 ± 0.27 | 1.00 | 1.62 | 2.68 | 1.00 |
