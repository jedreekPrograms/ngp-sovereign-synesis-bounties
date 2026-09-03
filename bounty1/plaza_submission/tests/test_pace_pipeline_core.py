"""Behavioral tests for the self-contained Plaza scoring core."""

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pace_pipeline_core import (  # noqa: E402
    HistoneMeasurement,
    PaceMeasurement,
    PairedObservation,
    PipelineEvidence,
    SampleKey,
    compute_primary_correlations,
    evidence_is_ready,
    pair_measurements,
    pearson,
    validate_evidence,
)


def key(condition, replicate):
    """Build a compact sample key for tests."""

    return SampleKey(condition, replicate)


def complete_rows():
    """Return four complete observations with independently supplied values."""

    pace = [
        PaceMeasurement(key("WT", 1), 0.88),
        PaceMeasurement(key("SIRT1", 1), 0.96),
        PaceMeasurement(key("SIRT2", 1), 1.07),
        PaceMeasurement(key("SIRT3", 1), 1.14),
    ]
    histone = []
    for sample_key, h3k9ac, h3k56ac in (
        (key("WT", 1), -0.2, -0.1),
        (key("SIRT1", 1), 0.1, 0.2),
        (key("SIRT2", 1), 0.7, 0.6),
        (key("SIRT3", 1), 1.0, 1.1),
    ):
        histone.extend(
            (
                HistoneMeasurement(sample_key, "H3K9ac", h3k9ac),
                HistoneMeasurement(sample_key, "H3K56ac", h3k56ac),
            )
        )
    return pace, histone


def ready_evidence():
    """Return non-biological pipeline evidence that satisfies quality gates."""

    return PipelineEvidence(
        paired_observations=8,
        mapq_threshold=30,
        peak_fdr=0.01,
        source_is_public=True,
        correlations_computed=True,
        synthetic_values_used=False,
    )


def test_pairing_requires_same_condition_and_replicate():
    pace, histone = complete_rows()
    paired = pair_measurements(pace, histone)
    assert len(paired) == 4
    assert paired[0].key == key("SIRT1", 1)
    assert paired[-1].key == key("WT", 1)


def test_pairing_excludes_incomplete_observation():
    pace, histone = complete_rows()
    histone = [row for row in histone if not (row.key == key("SIRT2", 1) and row.mark == "H3K56ac")]
    paired = pair_measurements(pace, histone)
    assert len(paired) == 3
    assert key("SIRT2", 1) not in {row.key for row in paired}


def test_pairing_ignores_unrequested_histone_mark():
    pace, histone = complete_rows()
    histone.append(HistoneMeasurement(key("WT", 1), "H3K27ac", 99.0))
    assert len(pair_measurements(pace, histone)) == 4


def test_duplicate_pace_measurement_is_rejected():
    pace, histone = complete_rows()
    pace.append(PaceMeasurement(key("WT", 1), 0.91))
    with pytest.raises(ValueError, match="duplicate DunedinPACE"):
        pair_measurements(pace, histone)


def test_duplicate_histone_measurement_is_rejected():
    pace, histone = complete_rows()
    histone.append(HistoneMeasurement(key("WT", 1), "H3K9ac", 0.3))
    with pytest.raises(ValueError, match="duplicate histone"):
        pair_measurements(pace, histone)


def test_nonfinite_pace_is_rejected():
    pace, histone = complete_rows()
    pace[0] = PaceMeasurement(key("WT", 1), math.nan)
    with pytest.raises(ValueError, match="non-finite DunedinPACE"):
        pair_measurements(pace, histone)


def test_nonfinite_histone_is_rejected():
    pace, histone = complete_rows()
    histone[0] = HistoneMeasurement(key("WT", 1), "H3K9ac", math.inf)
    with pytest.raises(ValueError, match="non-finite histone"):
        pair_measurements(pace, histone)


def test_pearson_perfect_positive_relationship():
    result = pearson([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 6.0, 8.0])
    assert result.n == 4
    assert result.pearson_r == pytest.approx(1.0)


def test_pearson_perfect_negative_relationship():
    result = pearson([1.0, 2.0, 3.0], [9.0, 6.0, 3.0])
    assert result.pearson_r == pytest.approx(-1.0)


def test_pearson_requires_equal_lengths():
    with pytest.raises(ValueError, match="equal length"):
        pearson([1.0, 2.0, 3.0], [1.0, 2.0])


def test_pearson_requires_three_observations():
    with pytest.raises(ValueError, match="at least three"):
        pearson([1.0, 2.0], [1.0, 2.0])


def test_pearson_rejects_zero_variance():
    with pytest.raises(ValueError, match="zero variance"):
        pearson([1.0, 1.0, 1.0], [2.0, 3.0, 4.0])


def test_primary_correlations_are_computed_for_both_marks():
    pace, histone = complete_rows()
    correlations = compute_primary_correlations(pair_measurements(pace, histone))
    assert set(correlations) == {"H3K9ac_vs_DunedinPACE", "H3K56ac_vs_DunedinPACE"}
    assert correlations["H3K9ac_vs_DunedinPACE"].n == 4
    assert correlations["H3K56ac_vs_DunedinPACE"].n == 4


def test_primary_correlations_require_three_complete_observations():
    rows = [
        PairedObservation(key("WT", 1), 0.9, 0.0, 0.0),
        PairedObservation(key("SIRT1", 1), 1.0, 0.2, 0.3),
    ]
    with pytest.raises(ValueError, match="three complete"):
        compute_primary_correlations(rows)


def test_ready_evidence_has_no_validation_problems():
    evidence = ready_evidence()
    assert validate_evidence(evidence) == ()
    assert evidence_is_ready(evidence) is True


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("paired_observations", 2, "insufficient paired"),
        ("mapq_threshold", 29, "MAPQ threshold"),
        ("peak_fdr", 0.05, "peak FDR"),
        ("source_is_public", False, "not publicly traceable"),
        ("correlations_computed", False, "not computed"),
        ("synthetic_values_used", True, "synthetic production"),
    ),
)
def test_evidence_rejects_each_failed_quality_gate(field, value, message):
    values = ready_evidence().__dict__.copy()
    values[field] = value
    problems = validate_evidence(PipelineEvidence(**values))
    assert any(message in problem for problem in problems)


def test_evidence_can_use_stricter_caller_thresholds():
    evidence = ready_evidence()
    problems = validate_evidence(
        evidence,
        minimum_pairs=6,
        minimum_mapq=30,
        maximum_peak_fdr=0.02,
    )
    assert problems == ()
