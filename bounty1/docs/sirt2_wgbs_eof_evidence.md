# SIRT2 WGBS checksum-verified EOF evidence

GitHub Actions run `33499862981` intentionally probed source-pair range
`450,000,000..499,999,999` for both SIRT2 WGBS biological replicates after
parts 1-9 had already covered the preceding non-overlapping 50M ranges.

The run is displayed by GitHub as `failure`, but the failure is **not** a
source-integrity or analysis failure.  For both samples the complete R1/R2
source stream finished, both archive-provided MD5 checks passed, and the
measured source pair count proved that `pair_start=450000000` was beyond EOF.
The range selector therefore emitted zero records and the old workflow treated
that valid EOF condition as an error.  No part10 candidate checkpoint exists
and no part10 retry is required.

## SIRT2_rep1_WGBS / HRR1202742

Measured by job `99830453480`:

- complete paired FASTQ records: **414,293,171**
- requested part10 start: `450,000,000`
- selected records: `0`
- R1 MD5 verified: `4acd4ee5c04819c9cafdba2e13f97e93`
- R2 MD5 verified: `3c8e437a40220f08ea729a28204ae505`
- failure-diagnostics artifact: `bounty1-wgbs-segmented-SIRT2_rep1_WGBS-part10-diagnostics`
- diagnostics artifact ZIP digest: `sha256:c0469580ff3e2746c7e0bc06af7d6e186db66c8e666f49d08ed7cc15fde8cd5b`

The durable part9 range (`400M..450M`) therefore contains the terminal source
segment `400,000,000..414,293,170`; no source records exist after it.

## SIRT2_rep2_WGBS / HRR1202743

Measured by job `99830453153`:

- complete paired FASTQ records: **410,082,682**
- requested part10 start: `450,000,000`
- selected records: `0`
- R1 MD5 verified: `bc3669802e92a739d37333e05f662bab`
- R2 MD5 verified: `bfe636eef625b7146169a4730aa350d6`
- failure-diagnostics artifact: `bounty1-wgbs-segmented-SIRT2_rep2_WGBS-part10-diagnostics`
- diagnostics artifact ZIP digest: `sha256:2d4f935a5b37c0af62974a5fe9a26d512cb1628d8a391d2b1daa51f175969e92`

The durable part9 range (`400M..450M`) therefore contains the terminal source
segment `400,000,000..410,082,681`; no source records exist after it.

## Workflow correction

`bin/screen_wgbs_pace_candidates.sh` now treats this specific condition as a
successful EOF checkpoint **only** when:

1. the complete source stream succeeds,
2. published R1 and R2 MD5 values both verify,
3. the measured pair-range log reports `selected=0`, and
4. measured `pair_start >= source_total_pairs`.

Other zero-output or failed-stream conditions remain fail-closed.
