# Unrelaxed rigid-rotor torsion MC (300 K)

| comparison | metric | median log(a/b) | a better / b better | p (Holm) |
|---|---|---|---|---|
| perlin vs gauss-corr | esjd | +0.033 | 33 / 21 | 0.446 (1) |
| perlin vs gauss-corr | transitions_per_1000 | -0.097 | 22 / 32 | 0.142 (0.71) |
| perlin vs gauss-corr | acceptance | -0.024 | 20 / 34 | 0.0272 (0.24) |
| perlin-gaussmarg vs perlin | esjd | -0.013 | 21 / 33 | 0.359 (1) |
| perlin-gaussmarg vs perlin | transitions_per_1000 | +0.084 | 32 / 20 | 0.0282 (0.24) |
| perlin-gaussmarg vs perlin | acceptance | +0.019 | 36 / 18 | 0.0501 (0.3) |
| gauss-perlinmarg vs gauss-corr | esjd | +0.025 | 31 / 23 | 0.359 (1) |
| gauss-perlinmarg vs gauss-corr | transitions_per_1000 | -0.022 | 24 / 28 | 0.168 (0.71) |
| gauss-perlinmarg vs gauss-corr | acceptance | -0.013 | 17 / 37 | 0.0388 (0.27) |
| gauss-corr vs gauss-iid | esjd | +0.073 | 41 / 13 | 1.57e-05 (0.00017) |
| gauss-corr vs gauss-iid | transitions_per_1000 | +0.196 | 41 / 12 | 7.72e-05 (0.00077) |
| gauss-corr vs gauss-iid | acceptance | +0.205 | 52 / 2 | 2.62e-09 (3.1e-08) |

| molecule | sigma | method | acceptance | ESJD (rad^2) | transitions/1000 | <E> |
|---|---|---|---|---|---|---|
| decane | 5 | gauss-iid | 0.543 | 0.0260 | 0.90 | -1.80 |
| decane | 5 | gauss-corr | 0.587 | 0.0261 | 1.89 | -1.64 |
| decane | 5 | perlin | 0.590 | 0.0274 | 1.63 | -1.29 |
| decane | 5 | perlin-gaussmarg | 0.586 | 0.0264 | 1.44 | -1.49 |
| decane | 5 | gauss-perlinmarg | 0.585 | 0.0268 | 1.77 | -1.68 |
| decane | 10 | gauss-iid | 0.246 | 0.0400 | 0.88 | -1.70 |
| decane | 10 | gauss-corr | 0.308 | 0.0432 | 1.33 | -1.42 |
| decane | 10 | perlin | 0.303 | 0.0450 | 1.03 | -1.38 |
| decane | 10 | perlin-gaussmarg | 0.308 | 0.0437 | 0.88 | -2.11 |
| decane | 10 | gauss-perlinmarg | 0.304 | 0.0447 | 1.37 | -1.37 |
| decane | 20 | gauss-iid | 0.043 | 0.0189 | 1.13 | -1.41 |
| decane | 20 | gauss-corr | 0.081 | 0.0275 | 0.92 | -1.54 |
| decane | 20 | perlin | 0.073 | 0.0263 | 0.88 | -1.35 |
| decane | 20 | perlin-gaussmarg | 0.078 | 0.0266 | 0.82 | -1.74 |
| decane | 20 | gauss-perlinmarg | 0.077 | 0.0273 | 0.85 | -1.51 |
| tetraglyme | 5 | gauss-iid | 0.496 | 0.0422 | 9.94 | 99.64 |
| tetraglyme | 5 | gauss-corr | 0.516 | 0.0410 | 13.17 | 100.95 |
| tetraglyme | 5 | perlin | 0.519 | 0.0426 | 14.38 | 100.75 |
| tetraglyme | 5 | perlin-gaussmarg | 0.526 | 0.0426 | 12.32 | 101.20 |
| tetraglyme | 5 | gauss-perlinmarg | 0.517 | 0.0418 | 15.17 | 100.89 |
| tetraglyme | 10 | gauss-iid | 0.188 | 0.0567 | 10.24 | 99.57 |
| tetraglyme | 10 | gauss-corr | 0.235 | 0.0612 | 11.36 | 100.13 |
| tetraglyme | 10 | perlin | 0.223 | 0.0615 | 10.63 | 99.96 |
| tetraglyme | 10 | perlin-gaussmarg | 0.228 | 0.0606 | 12.47 | 100.22 |
| tetraglyme | 10 | gauss-perlinmarg | 0.228 | 0.0618 | 10.07 | 99.51 |
| tetraglyme | 20 | gauss-iid | 0.015 | 0.0140 | 3.06 | 100.09 |
| tetraglyme | 20 | gauss-corr | 0.040 | 0.0275 | 4.78 | 100.06 |
| tetraglyme | 20 | perlin | 0.034 | 0.0260 | 4.45 | 100.48 |
| tetraglyme | 20 | perlin-gaussmarg | 0.037 | 0.0269 | 5.12 | 100.69 |
| tetraglyme | 20 | gauss-perlinmarg | 0.036 | 0.0262 | 4.84 | 100.23 |
| ala3 | 5 | gauss-iid | 0.495 | 0.0196 | 111.58 | 2.85 |
| ala3 | 5 | gauss-corr | 0.545 | 0.0206 | 106.69 | 3.46 |
| ala3 | 5 | perlin | 0.536 | 0.0208 | 98.86 | 3.72 |
| ala3 | 5 | perlin-gaussmarg | 0.534 | 0.0202 | 112.33 | 3.36 |
| ala3 | 5 | gauss-perlinmarg | 0.549 | 0.0213 | 95.20 | 3.65 |
| ala3 | 10 | gauss-iid | 0.214 | 0.0284 | 79.97 | 2.58 |
| ala3 | 10 | gauss-corr | 0.233 | 0.0268 | 89.44 | 2.64 |
| ala3 | 10 | perlin | 0.233 | 0.0284 | 88.92 | 2.45 |
| ala3 | 10 | perlin-gaussmarg | 0.241 | 0.0284 | 88.11 | 2.82 |
| ala3 | 10 | gauss-perlinmarg | 0.246 | 0.0302 | 87.09 | 3.13 |
| ala3 | 20 | gauss-iid | 0.041 | 0.0157 | 22.70 | 2.52 |
| ala3 | 20 | gauss-corr | 0.066 | 0.0210 | 33.77 | 2.96 |
| ala3 | 20 | perlin | 0.056 | 0.0188 | 29.77 | 2.86 |
| ala3 | 20 | perlin-gaussmarg | 0.061 | 0.0196 | 32.67 | 2.75 |
| ala3 | 20 | gauss-perlinmarg | 0.047 | 0.0149 | 26.02 | 2.53 |
