# Predicted vs. experimental conformations

16 molecules with experimental 3-D structures; 10 seeds x 300 MMFF94 relaxations per method; sigma = 60 deg, L = 4.
RMSD: heavy atoms, symmetry-aware, proper rotations only.

## Structures

| molecule | source | ID | resolution / R | torsions | heavy atoms | E(exp, H-opt) − E* | E(exp relaxed) − E* | RMSD exp→own minimum |
|---|---|---|---|---|---|---|---|---|
| imatinib | PDB | 1XBB | 1.57 Å | 7 | 37 | 60.3 | 5.8 | 0.93 |
| dasatinib | PDB | 6FNM | 1.16 Å | 7 | 33 | 55.6 | 14.2 | 0.55 |
| erlotinib | PDB | 9D5P | 1.50 Å | 10 | 29 | 23.5 | 8.8 | 0.37 |
| gefitinib | PDB | 9D5E | 1.50 Å | 8 | 31 | 40.4 | 7.0 | 0.31 |
| sildenafil | PDB | 1TBF | 1.30 Å | 7 | 33 | 41.0 | 1.4 | 0.74 |
| rosiglitazone | PDB | 9V8H | 1.39 Å | 7 | 25 | 16.1 | 4.5 | 0.60 |
| oseltamivir-carboxylate | PDB | 5NZ4 | 1.36 Å | 6 | 20 | 28.6 | 6.1 | 0.27 |
| ibuprofen | PDB | 3P6H | 1.15 Å | 4 | 15 | 25.5 | 0.6 | 0.73 |
| indomethacin | PDB | 7FCG | 1.19 Å | 5 | 25 | 28.8 | 0.8 | 0.67 |
| biotin | PDB | 2F01 | 0.85 Å | 5 | 16 | 42.8 | 9.7 | 0.49 |
| methotrexate | PDB | 3DAU | 1.50 Å | 9 | 33 | 89.0 | 11.6 | 0.74 |
| tolbutamide | COD | 7225537 | R = 0.027 | 5 | 18 | 13.1 | 2.7 | 0.79 |
| aspirin | COD | 7247819 | R = 0.052 | 3 | 13 | 8.2 | 1.9 | 0.08 |
| fenofibrate | COD | 2201994 | R = 0.052 | 7 | 25 | 8.1 | 0.2 | 0.55 |
| ketoprofen | COD | 2021652 | R = 0.039 | 4 | 19 | 21.4 | 0.9 | 0.31 |
| ibuprofen-crystal | COD | 2202426 | R = 0.038 | 4 | 15 | 18.6 | 0.9 | 0.21 |

## Methods (pooled over molecules and seeds)

| method | success ≤1.0 Å | success ≤0.5 Å | median min-RMSD | median min-RMSD (≤10 kcal) | median RMSD of lowest-E minimum | restricted mean steps to ≤1 Å |
|---|---|---|---|---|---|---|
| etkdg | 0.96 | 0.39 | 0.60 | 0.60 | 1.84 | 37 |
| reset | 0.90 | 0.35 | 0.66 | 0.66 | 2.02 | 62 |
| gauss-iid | 0.90 | 0.34 | 0.67 | 0.67 | 2.03 | 61 |
| gauss-corr | 0.90 | 0.29 | 0.72 | 0.73 | 1.87 | 61 |
| perlin | 0.91 | 0.34 | 0.67 | 0.67 | 1.94 | 66 |

## Paired tests over (molecule, seed)

| metric | comparison | mean diff [95% CI] | a better / b better | p (Holm) |
|---|---|---|---|---|
| min_rmsd | perlin − gauss-corr | -0.010 [-0.039, +0.017] | 74 / 77 | 0.837 (1) |
| min_rmsd | perlin − gauss-iid | +0.000 [-0.023, +0.023] | 77 / 74 | 0.717 (1) |
| min_rmsd | perlin − etkdg | +0.074 [+0.043, +0.106] | 75 / 76 | 0.077 (0.31) |
| min_rmsd | perlin − reset | +0.012 [-0.012, +0.037] | 75 / 76 | 0.389 (1) |
| rmsd_best | perlin − gauss-corr | +0.001 [-0.054, +0.054] | 86 / 69 | 0.34 (1) |
| rmsd_best | perlin − gauss-iid | -0.013 [-0.082, +0.059] | 85 / 70 | 0.486 (1) |
| rmsd_best | perlin − etkdg | +0.130 [+0.047, +0.220] | 60 / 96 | 0.00313 (0.013) |
| rmsd_best | perlin − reset | -0.043 [-0.109, +0.021] | 82 / 74 | 0.38 (1) |
| t_1_censored | perlin − gauss-corr | +5.487 [-8.925, +19.351] | 59 / 61 | 0.5 (1) |
| t_1_censored | perlin − gauss-iid | +5.419 [-7.170, +17.733] | 55 / 57 | 0.281 (0.84) |
| t_1_censored | perlin − etkdg | +28.962 [+13.618, +45.275] | 39 / 80 | 2.94e-05 (0.00012) |
| t_1_censored | perlin − reset | +3.956 [-7.576, +16.044] | 59 / 59 | 0.831 (1) |
| success_1 | perlin − gauss-corr | +0.006 [-0.044, +0.062] | 10 / 9 | 0.819 (1) |
| success_1 | perlin − gauss-iid | +0.006 [-0.031, +0.050] | 6 / 5 | 0.763 (1) |
| success_1 | perlin − etkdg | -0.056 [-0.106, -0.013] | 3 / 12 | 0.0201 (0.081) |
| success_1 | perlin − reset | +0.006 [-0.044, +0.056] | 9 / 8 | 0.808 (1) |

## Per molecule: success ≤1.0 Å / median min-RMSD

| molecule | etkdg | reset | gauss-iid | gauss-corr | perlin |
|---|---|---|---|---|---|
| imatinib | 0.8 / 0.87 | 1.0 / 0.78 | 0.9 / 0.85 | 0.7 / 0.93 | 0.9 / 0.95 |
| dasatinib | 1.0 / 0.66 | 0.3 / 1.17 | 0.1 / 1.30 | 0.5 / 1.01 | 0.2 / 1.21 |
| erlotinib | 1.0 / 0.64 | 1.0 / 0.67 | 1.0 / 0.71 | 1.0 / 0.67 | 1.0 / 0.64 |
| gefitinib | 1.0 / 0.38 | 0.9 / 0.77 | 1.0 / 0.77 | 0.9 / 0.83 | 1.0 / 0.80 |
| sildenafil | 1.0 / 0.78 | 1.0 / 0.74 | 1.0 / 0.73 | 1.0 / 0.73 | 1.0 / 0.73 |
| rosiglitazone | 1.0 / 0.60 | 0.9 / 0.60 | 0.7 / 0.75 | 0.8 / 0.89 | 0.6 / 0.88 |
| oseltamivir-carboxylate | 1.0 / 0.27 | 1.0 / 0.27 | 1.0 / 0.27 | 0.9 / 0.50 | 1.0 / 0.28 |
| ibuprofen | 1.0 / 0.73 | 1.0 / 0.73 | 1.0 / 0.73 | 1.0 / 0.73 | 1.0 / 0.73 |
| indomethacin | 1.0 / 0.61 | 1.0 / 0.61 | 1.0 / 0.61 | 1.0 / 0.61 | 1.0 / 0.61 |
| biotin | 1.0 / 0.45 | 1.0 / 0.45 | 1.0 / 0.45 | 1.0 / 0.62 | 1.0 / 0.45 |
| methotrexate | 0.6 / 0.99 | 0.3 / 1.10 | 0.7 / 0.95 | 0.6 / 0.86 | 0.8 / 0.85 |
| tolbutamide | 1.0 / 0.74 | 1.0 / 0.74 | 1.0 / 0.74 | 1.0 / 0.74 | 1.0 / 0.74 |
| aspirin | 1.0 / 0.08 | 1.0 / 0.08 | 1.0 / 0.08 | 1.0 / 0.08 | 1.0 / 0.08 |
| fenofibrate | 1.0 / 0.55 | 1.0 / 0.45 | 1.0 / 0.45 | 1.0 / 0.45 | 1.0 / 0.45 |
| ketoprofen | 1.0 / 0.31 | 1.0 / 0.31 | 1.0 / 0.31 | 1.0 / 0.31 | 1.0 / 0.31 |
| ibuprofen-crystal | 1.0 / 0.21 | 1.0 / 0.21 | 1.0 / 0.21 | 1.0 / 0.21 | 1.0 / 0.21 |
