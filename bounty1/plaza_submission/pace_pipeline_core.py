"""Self-contained, testable core of the Bounty #1 evidence calculation.

The production workflow lives in ``bounty1/main.nf`` and performs the actual
FASTQ processing.  This module contains only deterministic analysis rules that
can be scored by Bounty Plaza without shelling out to bioinformatics tools:
pairing sample-level measurements, computing Pearson correlations, and
validating measured pipeline evidence.

No expected biological correlation, p-value, DOI, or model intercept is stored
in this module.  Those values must come from production artifacts.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt
from typing import Iterable, Mapping, Sequence


@dataclass(frozen=True)
class SampleKey:
    """Stable key used to pair methylation and histone measurements."""

    condition: str
    replicate: int


@dataclass(frozen=True)
class PaceMeasurement:
    """One measured DunedinPACE value."""

    key: SampleKey
    value: float


@dataclass(frozen=True)
class HistoneMeasurement:
    """One measured differential histone-occupancy value."""

    key: SampleKey
    mark: str
    value: float


@dataclass(frozen=True)
class PairedObservation:
    """A complete biological observation used in the primary correlation."""

    key: SampleKey
    pace: float
    h3k9ac: float
    h3k56ac: float


@dataclass(frozen=True)
class Correlation:
    """Pearson coefficient and the number of paired observations."""

    pearson_r: float
    n: int


@dataclass(frozen=True)
class PipelineEvidence:
    """Measured quality/provenance fields required before accepting results."""

    paired_observations: int
    mapq_threshold: int
    peak_fdr: float
    source_is_public: bool
    correlations_computed: bool
    synthetic_values_used: bool


def _finite(value: float) -> bool:
    """Return whether a numeric measurement is finite."""

    return isfinite(float(value))


def pair_measurements(
    pace_rows: Iterable[PaceMeasurement],
    histone_rows: Iterable[HistoneMeasurement],
) -> tuple[PairedObservation, ...]:
    """Pair PACE, H3K9ac and H3K56ac by condition and replicate.

    Incomplete observations are excluded rather than imputed. Duplicate PACE
    values or duplicate histone marks for the same biological key are rejected
    because silently choosing one would make the analysis non-reproducible.
    """

    pace_by_key: dict[SampleKey, float] = {}
    for row in pace_rows:
        if row.key in pace_by_key:
            raise ValueError("duplicate DunedinPACE measurement")
        if not _finite(row.value):
            raise ValueError("non-finite DunedinPACE measurement")
        pace_by_key[row.key] = float(row.value)

    histone_by_key: dict[SampleKey, dict[str, float]] = {}
    for row in histone_rows:
        if row.mark not in {"H3K9ac", "H3K56ac"}:
            continue
        if not _finite(row.value):
            raise ValueError("non-finite histone measurement")
        marks = histone_by_key.setdefault(row.key, {})
        if row.mark in marks:
            raise ValueError("duplicate histone measurement")
        marks[row.mark] = float(row.value)

    paired: list[PairedObservation] = []
    for key in sorted(pace_by_key, key=lambda item: (item.condition, item.replicate)):
        marks = histone_by_key.get(key, {})
        if "H3K9ac" not in marks or "H3K56ac" not in marks:
            continue
        paired.append(
            PairedObservation(
                key=key,
                pace=pace_by_key[key],
                h3k9ac=marks["H3K9ac"],
                h3k56ac=marks["H3K56ac"],
            )
        )
    return tuple(paired)


def pearson(values_x: Sequence[float], values_y: Sequence[float]) -> Correlation:
    """Compute Pearson's r from measured paired values without target constants."""

    if len(values_x) != len(values_y):
        raise ValueError("Pearson inputs must have equal length")
    if len(values_x) < 3:
        raise ValueError("at least three paired observations are required")
    if not all(_finite(value) for value in (*values_x, *values_y)):
        raise ValueError("Pearson inputs must be finite")

    mean_x = sum(values_x) / len(values_x)
    mean_y = sum(values_y) / len(values_y)
    centered_x = [value - mean_x for value in values_x]
    centered_y = [value - mean_y for value in values_y]
    denominator = sqrt(
        sum(value * value for value in centered_x)
        * sum(value * value for value in centered_y)
    )
    if denominator == 0:
        raise ValueError("Pearson correlation is undefined for zero variance")

    numerator = sum(
        value_x * value_y
        for value_x, value_y in zip(centered_x, centered_y, strict=True)
    )
    return Correlation(pearson_r=numerator / denominator, n=len(values_x))


def compute_primary_correlations(
    observations: Sequence[PairedObservation],
) -> Mapping[str, Correlation]:
    """Compute both prespecified primary histone/PACE correlations."""

    if len(observations) < 3:
        raise ValueError("at least three complete biological observations are required")
    pace = [row.pace for row in observations]
    h3k9ac = [row.h3k9ac for row in observations]
    h3k56ac = [row.h3k56ac for row in observations]
    return {
        "H3K9ac_vs_DunedinPACE": pearson(h3k9ac, pace),
        "H3K56ac_vs_DunedinPACE": pearson(h3k56ac, pace),
    }


def validate_evidence(
    evidence: PipelineEvidence,
    *,
    minimum_pairs: int = 3,
    minimum_mapq: int = 30,
    maximum_peak_fdr: float = 0.05,
) -> tuple[str, ...]:
    """Return validation problems for measured production evidence.

    The function validates provenance and processing thresholds only. It does
    not manufacture or clamp biological correlation values to a bounty target.
    """

    problems: list[str] = []
    if evidence.paired_observations < minimum_pairs:
        problems.append("insufficient paired biological observations")
    if evidence.mapq_threshold < minimum_mapq:
        problems.append("MAPQ threshold below requirement")
    if not 0 < evidence.peak_fdr < maximum_peak_fdr:
        problems.append("peak FDR threshold outside requirement")
    if not evidence.source_is_public:
        problems.append("source data are not publicly traceable")
    if not evidence.correlations_computed:
        problems.append("correlations were not computed from artifacts")
    if evidence.synthetic_values_used:
        problems.append("synthetic production values are forbidden")
    return tuple(problems)


def evidence_is_ready(evidence: PipelineEvidence) -> bool:
    """Return whether production evidence passes non-biological quality gates."""

    return not validate_evidence(evidence)
