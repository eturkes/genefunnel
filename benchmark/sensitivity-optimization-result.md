<!-- Assisted-by: OpenAI Codex. -->

# Exact sorted-prefix adoption result

**Decision:** adopt the unexported exact sorted-prefix path; retain the brute
implementation as the executable oracle and fallback.

Candidate `e5c013fcb0f10e2eb4f0431c1cbd09cd40f73954` was installed from a
clean Git archive in an isolated library. The archive SHA-256 is
`11b7892eb9d60bea77ee63301673da0bae5e0a2744718952a88301f6ede79224`;
the installed-package manifest MD5 is `7f30f8246b98f6a1262a55b2f5217c09`.

The fixed E-P-1.0.0 fixture retained MD5
`32fa8d4843d024765a4adfc676793dbe`. Its candidate output MD5 was
`3d9635e779a9ed1eee453a2a04596369`, exactly the previously committed brute
output. Package tests additionally require exact object identity - including
signed limbs, denominator, total, exponent, ties, and statuses - over fixed,
wide-exponent, randomized, and size-128 cells, plus the existing storage and
backend suite.

The isolated call elapsed 15.128 seconds, recorded only as exploratory context.
It is not a replicated performance result and sets `performance_claim = FALSE`.
No optimization result establishes predictive reliability or supports a public
API. The machine-readable companion is
[`sensitivity-optimization-result.tsv`](sensitivity-optimization-result.tsv).
