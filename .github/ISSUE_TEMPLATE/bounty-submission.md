---
name: Bounty Submission
about: Submit your work for one of the NGP Sovereign Synesis bounties
title: "[SUBMISSION] Bounty #X — Your Name / Team"
labels: submission, under-review
assignees: mister3ai-cmyk
---

## Bounty Reference

- **Bounty #**: <!-- 1, 2, or 3 -->
- **Prize**: <!-- $15,000 / $25,000 / $20,000 USDC -->
- **Submitter(s)**: <!-- GitHub handle(s) -->
- **Ethereum address for payment**: `0x...`

## Repository / Artifact Links

- Code repository: <!-- URL -->
- Data deposit DOI: <!-- Zenodo / GEO DOI -->
- Preprint (if applicable): <!-- arXiv URL -->
- Video demo (if applicable): <!-- YouTube / Loom URL -->

## CI Results

Paste your local `pytest` run output here:

```
pytest tests/test_bountyX_*.py -v
```

```
[paste output here]
```

## Checklist

- [ ] All CI tests pass locally
- [ ] `results/manifest.json` (or equivalent) is included in the repository
- [ ] Docker image / Docker Compose stack is publicly accessible or image is attached
- [ ] Report / specification document is attached or linked
- [ ] Data deposited with persistent DOI
- [ ] I confirm this is original work and I hold the rights to license it under CC BY 4.0

## Notes for Reviewers

<!-- Any additional context, known limitations, or requests for partial credit -->
