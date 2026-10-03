# ANI-2x (GPU) as a second judge

## Conformational isomers: MMFF94 vs ANI-2x vs gas-phase experiment

| molecule | quantity | MMFF94 | ANI-2x | experiment |
|---|---|---|---|---|
| n-butane | gauche - anti | +0.78 | +0.55 | +0.66 ± 0.02 |
| 1,2-difluoroethane | gauche - anti | -0.63 | -0.10 | -0.80 ± 0.50 |
| 1,2-dichloroethane | gauche - anti | +1.24 | +1.14 | +1.05 ± 0.10 |
| n-pentane | tg - tt | +0.83 | +0.50 | +0.62 ± 0.01 |
| n-pentane | g+g+ - tt | +1.42 | +0.98 | +0.94 ± 0.02 |

MAE vs experiment: MMFF94 0.23, ANI-2x 0.21 kcal/mol (same molecules)

## Re-ranking the discovered low-energy minima

| molecule | minima | Spearman(MMFF, ANI) | rank of MMFF global min under ANI | MMFF dE of ANI global min | median geometry shift (A) |
|---|---|---|---|---|---|
| ala3 | 9 | 0.23 | 3 | 2.27 | 0.16 |
| decane | 193 | 0.90 | 0 | 0.00 | 0.02 |
| imatinib | 88 | 0.16 | 18 | 1.58 | 0.15 |
| tetraglyme | 56 | 0.59 | 0 | 0.00 | 0.08 |
