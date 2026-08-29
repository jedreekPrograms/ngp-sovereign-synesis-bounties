# Production launch controls

Heavy Bounty #1 production workflows are decoupled from pull-request synchronize
events.

For segmented WGBS, first verify that no overlapping shard is active and that all
earlier ranges have durable, checksum-verified artifacts. Then create or change
`wgbs-launch-token.txt` in a dedicated commit. That path is the only automatic
push trigger for the armed WGBS production matrix.

Do not change the token merely to re-run a failed shard. Diagnose the failure and
perform the smallest recovery that preserves existing checkpoints.
