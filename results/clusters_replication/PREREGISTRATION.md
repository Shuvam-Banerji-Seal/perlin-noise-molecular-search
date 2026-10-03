# Pre-registered replication + mechanism test (written 2026-10-03T01:48Z, before any run)

History. LJ38 basin hopping (3000 relaxations, s = 0.3 sigma, kT = 0.8):
- exploratory, seeds 0-29: perlin@2 17/30 vs gauss-field@2 11/30 (Fisher p = 0.195);
- pre-registered confirmation, seeds 30-129: 45/100 vs 35/100 (Fisher p = 0.194) -> NOT confirmed.
The direction was the same both times. The confirmation had about 30% power for a 10-point effect;
about 380 runs per arm are needed for 80% power. At L = 2 sigma the Perlin displacement field is almost
uncorrelated between atoms (mean inter-atomic correlation 0.025), so any effect must come from its
light-tailed marginal distribution (excess kurtosis about -0.32).

H3 (primary, replication): P_hit(perlin@2) != P_hit(gauss-field@2).
  Design: 300 fresh seeds (130-429), same protocol. Two-sided Fisher exact test, alpha = 0.05.
H4 (secondary, mechanism): an i.i.d. kick drawn from Perlin's marginal distribution (iid-perlinmarg,
  no spatial correlation) on the same 300 seeds. Two-sided Fisher tests vs gauss-field@2 and vs
  perlin@2, reported with Holm correction over these two.
Descriptive only (not a test): pooled seeds 30-429 (n = 400) for perlin@2 vs gauss-field@2.
Decision: H3 is "replicated" only if p < 0.05 with Perlin ahead; otherwise "not replicated".
