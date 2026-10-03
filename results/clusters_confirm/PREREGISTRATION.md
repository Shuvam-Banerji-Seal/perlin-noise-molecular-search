# Pre-registered confirmatory test (written 2026-10-02T21:19Z, before any confirmatory run)

Exploratory finding (results/clusters, seeds 0-29): on LJ38, basin hopping with 3000 relaxations
(s = 0.3 sigma, kT = 0.8) found the known global minimum in 17/30 runs with perlin@2 vs 11/30 with its
matched Gaussian field gauss-field@2 (Fisher p = 0.195; not significant). Direction also seen at L = 4 (10/30 vs 6/30).
This was one favourable result among many comparisons, so it is a hypothesis, not a finding.

Primary hypothesis H1: P_hit(perlin@2) != P_hit(gauss-field@2) on LJ38 under the same protocol.
Design: 100 fresh seeds (30-129, disjoint from the exploratory seeds), same budget/step/kT, same start per seed.
Primary test: two-sided Fisher exact test on success counts, alpha = 0.05.
Secondary (descriptive only): paired difference in final gap and in censored time-to-GM.
Decision rule: if p >= 0.05, the LJ38 hint is reported as not confirmed.
