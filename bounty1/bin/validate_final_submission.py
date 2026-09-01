#!/usr/bin/env python3
"""Strict final-submission validator for Bounty #1.

This validator is intentionally stricter than the upstream pytest harness.  It
checks the literal scientific acceptance criteria as well as the Definition of
Done evidence that can be verified locally.  It never substitutes acceptance
targets for measured values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

REFERENCE_INTERCEPT = 51.024577
INTERCEPT_TOLERANCE = 0.001
MIN_PEARSON_R = 0.92
MAX_P_VALUE = 0.01
MAX_FDR = 0.05
MIN_MAPQ = 30
MIN_REPORT_PAGES = 8


def load_json(path: Path) -> dict:
    if not path.is_file():
        raise ValueError(f"missing required file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing required file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(value, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} is not finite")
    return number


def md5sum(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdf_pages(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing/empty supplementary report: {path}")
    try:
        from pypdf import PdfReader
    except ImportError as exc:  # pragma: no cover - dependency installed in final CI
        raise ValueError("pypdf is required to validate report page count") from exc
    return len(PdfReader(str(path)).pages)


def validate(
    manifest_path: Path,
    correlations_path: Path,
    scores_path: Path,
    occupancy_path: Path,
    pace_qc_path: Path,
    report_path: Path,
    min_paired: int,
) -> dict:
    manifest = load_json(manifest_path)
    correlations = load_json(correlations_path)
    scores = load_csv(scores_path)
    occupancy = load_csv(occupancy_path)
    pace_qc = load_csv(pace_qc_path)

    errors: list[str] = []

    # Upstream acceptance-test contract.
    try:
        intercept = finite(manifest["dunedinpace"]["intercept"], "DunedinPACE intercept")
        if abs(intercept - REFERENCE_INTERCEPT) > INTERCEPT_TOLERANCE:
            errors.append(
                f"DunedinPACE intercept {intercept:.6f} does not meet "
                f"{REFERENCE_INTERCEPT:.6f} ± {INTERCEPT_TOLERANCE}"
            )
    except (KeyError, TypeError, ValueError) as exc:
        errors.append(f"invalid DunedinPACE intercept: {exc}")

    n_paired = int(correlations.get("n_paired", 0))
    if n_paired < min_paired:
        errors.append(f"only {n_paired} fully paired observations; require >= {min_paired}")

    for mark in ("H3K9ac", "H3K56ac"):
        key = f"{mark}_vs_DunedinPACE"
        try:
            result = correlations[key]
            r_value = finite(result["pearson_r"], f"{mark} Pearson r")
            p_value = finite(result["p_value"], f"{mark} p-value")
            n_value = int(result.get("n", n_paired))
            if r_value <= MIN_PEARSON_R:
                errors.append(f"{mark} Pearson r={r_value:.6f} <= {MIN_PEARSON_R}")
            if p_value >= MAX_P_VALUE:
                errors.append(f"{mark} p={p_value:.6g} >= {MAX_P_VALUE}")
            if n_value < min_paired:
                errors.append(f"{mark} correlation has n={n_value}; require >= {min_paired}")
        except (KeyError, TypeError, ValueError) as exc:
            errors.append(f"invalid {mark} correlation: {exc}")

    try:
        fdr = finite(manifest["peak_calling"]["fdr"], "peak FDR")
        if fdr >= MAX_FDR:
            errors.append(f"peak FDR={fdr} >= {MAX_FDR}")
    except (KeyError, TypeError, ValueError) as exc:
        errors.append(f"invalid peak FDR: {exc}")

    try:
        mapq = int(manifest["alignment"]["mapq_threshold"])
        if mapq < MIN_MAPQ:
            errors.append(f"MAPQ threshold={mapq} < {MIN_MAPQ}")
    except (KeyError, TypeError, ValueError) as exc:
        errors.append(f"invalid MAPQ threshold: {exc}")

    doi = str(manifest.get("data_deposit_doi", "")).strip()
    if not doi.startswith("10.") or "/" not in doi:
        errors.append(f"missing/invalid persistent DOI: {doi!r}")

    provenance = manifest.get("provenance", {})
    if provenance.get("synthetic_values_used") is not False:
        errors.append("manifest must explicitly state provenance.synthetic_values_used=false")
    if provenance.get("correlations_computed") is not True:
        errors.append("manifest must explicitly state provenance.correlations_computed=true")

    # Docker/Singularity checksum criterion: require a real exported image file,
    # not merely a digest string detached from an artifact.
    image_path_raw = str(manifest.get("docker_image_path", "")).strip()
    expected_md5 = str(manifest.get("docker_image_md5", "")).strip().lower()
    if not image_path_raw or not expected_md5:
        errors.append("docker_image_path/docker_image_md5 missing from final manifest")
    else:
        image_path = Path(image_path_raw)
        if not image_path.is_file():
            errors.append(f"exported Docker image not found: {image_path}")
        else:
            actual_md5 = md5sum(image_path)
            if actual_md5 != expected_md5:
                errors.append(
                    f"Docker image MD5 mismatch: manifest={expected_md5}, actual={actual_md5}"
                )

    pages = 0
    try:
        pages = pdf_pages(report_path)
        if pages < MIN_REPORT_PAGES:
            errors.append(f"supplementary report has {pages} pages; require >= {MIN_REPORT_PAGES}")
    except ValueError as exc:
        errors.append(str(exc))

    if len(scores) < min_paired:
        errors.append(f"DunedinPACE score table has {len(scores)} rows; require >= {min_paired}")
    if len(pace_qc) < min_paired:
        errors.append(f"PACE QC table has {len(pace_qc)} rows; require >= {min_paired}")

    # Each paired observation needs both histone marks.  WT rows may additionally
    # be present solely as the prespecified centering baseline.
    keys_by_mark: dict[str, set[tuple[str, int]]] = {"H3K9ac": set(), "H3K56ac": set()}
    for row in occupancy:
        mark = row.get("mark", "")
        if mark not in keys_by_mark:
            continue
        try:
            keys_by_mark[mark].add((row["condition"], int(row["replicate"])))
        except (KeyError, ValueError):
            errors.append(f"invalid occupancy metadata row: {row}")
    common_occ = keys_by_mark["H3K9ac"] & keys_by_mark["H3K56ac"]
    if len(common_occ) < min_paired:
        errors.append(
            f"only {len(common_occ)} condition/replicate keys have both histone marks; "
            f"require >= {min_paired}"
        )
    for mark in ("H3K9ac", "H3K56ac"):
        if not any(condition == "WT" for condition, _ in keys_by_mark[mark]):
            errors.append(f"missing WT baseline occupancy for {mark}")

    summary = {
        "ok": not errors,
        "errors": errors,
        "n_paired": n_paired,
        "scores": len(scores),
        "pace_qc_rows": len(pace_qc),
        "report_pages": pages,
        "doi": doi,
    }
    if errors:
        raise ValueError("FINAL ACCEPTANCE NOT MET:\n- " + "\n- ".join(errors))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--correlations", required=True, type=Path)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--occupancy", required=True, type=Path)
    parser.add_argument("--pace-qc", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--min-paired", type=int, default=6)
    args = parser.parse_args()
    summary = validate(
        args.manifest,
        args.correlations,
        args.scores,
        args.occupancy,
        args.pace_qc,
        args.report,
        args.min_paired,
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
