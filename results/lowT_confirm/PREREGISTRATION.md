# Pre-registered confirmatory test 2 (written 2026-10-03T01:34Z, before any confirmatory run)

Exploratory finding (results/hunt/kT0.3, seeds 0-9): greedy basin hopping at kT = 0.3 kcal/mol,
sigma = 60 deg, 500 relaxations, on decane, tetraglyme and Ac-(Ala)3-NHMe. perlin(L=4) recovered more of
the grid's low-energy minima than its matched Gaussian gauss-corr(L=4): paired mean difference +0.047
[95% CI +0.013, +0.084], Wilcoxon p = 0.024, Holm-adjusted p = 0.07 within its family. The same contrast at
L = 8 was null (+0.004), and at kT = 3 null (+0.003). Many other tests were run, so this is a hypothesis.

Primary hypothesis H2: final recall of the grid low-energy set differs between perlin(L=4) and gauss-corr(L=4)
under the same protocol.
Design: 30 fresh seeds (10-39, disjoint from the exploratory seeds 0-9), same three molecules, same protocol;
paired by (molecule, seed). Primary test: two-sided Wilcoxon signed-rank on the 90 paired differences, alpha = 0.05.
Secondary (descriptive): paired bootstrap 95% CI; per-molecule means; time to the known global minimum.
Decision rule: p >= 0.05 -> reported as not confirmed.
