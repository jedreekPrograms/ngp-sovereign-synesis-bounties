# Disabled heavy GitHub Actions workflows

These workflows are retained byte-for-byte for provenance, but intentionally live
outside `.github/workflows/` so GitHub Actions cannot auto-run them.

Why: this long-lived draft PR changes many Bounty #1 paths. GitHub evaluates
`pull_request.paths` against the cumulative PR diff, so every later synchronize
event used to re-trigger historical WGBS benchmarks, smoke tests, completed
50M WGBS ranges, and the old WT chromatin workflow (including its matched INPUT).

Production replacements use explicit push control tokens or manual dispatch.
Restore an archived workflow only for a deliberate new validation experiment.
Never restore old production matrices without first checking durable artifacts
and active ranges.
