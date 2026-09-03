#!/usr/bin/env python3
"""Generate the Bounty #1 supplementary PDF strictly from measured outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.enums import TA_CENTER  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.platypus import (  # noqa: E402
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


DATASET_DOI = "10.1016/j.devcel.2024.02.008"
DATASET_ACCESSIONS = "HRA003336 / PRJCA012536; SIRT6 CUT&RUN: HRA005392"
PACE_SOURCE = "danbelsky/DunedinPACE@4b569983543e51d1022aecec9a25e694bb3a336a"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require_rows(rows: list[dict[str, str]], label: str, minimum: int = 1) -> None:
    if len(rows) < minimum:
        raise ValueError(f"{label} has {len(rows)} rows; need at least {minimum}")


def finite_float(value: str | float | int, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label} is not finite")
    return number


def plot_scores(scores: list[dict[str, str]], output: Path) -> None:
    labels = [row["sample_id"] for row in scores]
    values = [finite_float(row["dunedinpace"], "DunedinPACE") for row in scores]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(range(len(values)), values)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("DunedinPACE")
    ax.set_title("Measured DunedinPACE scores")
    fig.tight_layout()
    fig.savefig(output, dpi=170)
    plt.close(fig)


def plot_probe_qc(qc: list[dict[str, str]], output: Path) -> None:
    labels = [row["sample_id"] for row in qc]
    values = [finite_float(row["fraction"], "probe fraction") for row in qc]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.bar(range(len(values)), values)
    ax.axhline(0.8, linestyle="--", linewidth=1)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Fraction of required probes observed")
    ax.set_title("DunedinPACE probe coverage QC")
    fig.tight_layout()
    fig.savefig(output, dpi=170)
    plt.close(fig)


def paired_points(
    scores: list[dict[str, str]], occupancy: list[dict[str, str]], mark: str
) -> tuple[list[str], list[float], list[float]]:
    pace = {
        (row["condition"], int(row["replicate"])): finite_float(
            row["dunedinpace"], "DunedinPACE"
        )
        for row in scores
    }
    points = []
    for row in occupancy:
        if row["mark"] != mark:
            continue
        key = (row["condition"], int(row["replicate"]))
        if key in pace:
            points.append(
                (
                    f"{key[0]} r{key[1]}",
                    finite_float(row["differential_occupancy"], "occupancy"),
                    pace[key],
                )
            )
    if len(points) < 3:
        raise ValueError(f"Need at least 3 paired observations for {mark}")
    labels, x, y = zip(*points)
    return list(labels), list(x), list(y)


def plot_correlation(
    scores: list[dict[str, str]],
    occupancy: list[dict[str, str]],
    correlations: dict,
    mark: str,
    output: Path,
) -> None:
    labels, x, y = paired_points(scores, occupancy, mark)
    key = f"{mark}_vs_DunedinPACE"
    result = correlations[key]
    r_value = finite_float(result["pearson_r"], f"{mark} r")
    p_value = finite_float(result["p_value"], f"{mark} p")

    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    ax.scatter(x, y)
    for label, x_value, y_value in zip(labels, x, y):
        ax.annotate(label, (x_value, y_value), fontsize=7, xytext=(3, 3), textcoords="offset points")
    if len(x) >= 2 and max(x) != min(x):
        slope, intercept = _linear_fit(x, y)
        x_line = [min(x), max(x)]
        y_line = [slope * value + intercept for value in x_line]
        ax.plot(x_line, y_line, linewidth=1)
    ax.set_xlabel(f"Differential {mark} occupancy")
    ax.set_ylabel("DunedinPACE")
    ax.set_title(f"{mark} vs DunedinPACE: r={r_value:.4f}, p={p_value:.3g}")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def _linear_fit(x: list[float], y: list[float]) -> tuple[float, float]:
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0:
        return 0.0, y_mean
    slope = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y)) / denominator
    return slope, y_mean - slope * x_mean


def table(data: list[list[str]], widths: list[float] | None = None) -> Table:
    result = Table(data, colWidths=widths, repeatRows=1)
    result.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.grey),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return result


def add_page_number(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.drawRightString(200 * mm, 9 * mm, f"Page {doc.page}")
    canvas.restoreState()


def build_report(
    scores_path: Path,
    occupancy_path: Path,
    correlations_path: Path,
    manifest_path: Path,
    pace_qc_path: Path,
    output: Path,
) -> None:
    scores = read_csv(scores_path)
    occupancy = read_csv(occupancy_path)
    pace_qc = read_csv(pace_qc_path)
    correlations = read_json(correlations_path)
    manifest = read_json(manifest_path)
    require_rows(scores, "scores", 3)
    require_rows(occupancy, "occupancy", 6)
    require_rows(pace_qc, "PACE QC", 3)

    for mark in ("H3K9ac", "H3K56ac"):
        key = f"{mark}_vs_DunedinPACE"
        if key not in correlations:
            raise ValueError(f"Missing measured correlation: {key}")
        finite_float(correlations[key]["pearson_r"], f"{mark} r")
        finite_float(correlations[key]["p_value"], f"{mark} p")

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="TitleCenter",
            parent=styles["Title"],
            alignment=TA_CENTER,
            spaceAfter=14,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Small",
            parent=styles["BodyText"],
            fontSize=8.5,
            leading=11,
        )
    )

    doc = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        rightMargin=18 * mm,
        leftMargin=18 * mm,
        topMargin=18 * mm,
        bottomMargin=16 * mm,
        title="Bounty #1 - ChIP-seq and DunedinPACE supplementary report",
        author="jedreekPrograms",
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        scores_png = temp / "pace_scores.png"
        qc_png = temp / "pace_probe_qc.png"
        h3k9_png = temp / "h3k9_corr.png"
        h3k56_png = temp / "h3k56_corr.png"
        plot_scores(scores, scores_png)
        plot_probe_qc(pace_qc, qc_png)
        plot_correlation(scores, occupancy, correlations, "H3K9ac", h3k9_png)
        plot_correlation(scores, occupancy, correlations, "H3K56ac", h3k56_png)

        story = []
        story.extend(
            [
                Spacer(1, 35 * mm),
                Paragraph("Bounty #1", styles["TitleCenter"]),
                Paragraph(
                    "ChIP-seq & Methylation PACE Pipeline - Supplementary Evidence Report",
                    styles["TitleCenter"],
                ),
                Spacer(1, 12 * mm),
                Paragraph(
                    "This document is generated from measured pipeline outputs. Correlations, "
                    "p-values, probe coverage and manifest values are not substituted with "
                    "acceptance targets.",
                    styles["BodyText"],
                ),
                Spacer(1, 8 * mm),
                table(
                    [
                        ["Item", "Value"],
                        ["Primary dataset", DATASET_ACCESSIONS],
                        ["Source publication DOI", DATASET_DOI],
                        ["DunedinPACE implementation", PACE_SOURCE],
                        ["Samples with PACE scores", str(len(scores))],
                        ["MAPQ threshold", str(manifest.get("alignment", {}).get("mapq_threshold", "n/a"))],
                    ],
                    [48 * mm, 120 * mm],
                ),
                PageBreak(),
            ]
        )

        story.extend(
            [
                Paragraph("1. Data provenance and experimental design", styles["Heading1"]),
                Paragraph(
                    "The pipeline uses public human mesenchymal progenitor/stem-cell data from "
                    "HRA003336 (PRJCA012536) and SIRT6 CUT&RUN data from HRA005392. H3K9ac, "
                    "H3K56ac and WGBS observations are paired by condition and replicate. Matched "
                    "INPUT libraries are used for ChIP peak calling and matched IgG libraries are "
                    "used for the SIRT6 CUT&RUN locus definition.",
                    styles["BodyText"],
                ),
                Spacer(1, 5 * mm),
                Paragraph(
                    "Important scope limitation: DunedinPACE was developed and validated primarily "
                    "for blood DNA methylation. Application to hMPC WGBS is an off-domain analysis "
                    "and is interpreted as an experimental association, not a validated clinical "
                    "pace-of-aging estimate.",
                    styles["BodyText"],
                ),
                Spacer(1, 5 * mm),
                Paragraph(
                    "The source publication describes WGBS processing with fastp, hg19 mapping and "
                    "CpG depth >5. This workflow preserves the hg19 coordinate system and minimum "
                    "depth criterion while providing a reproducible containerized implementation.",
                    styles["BodyText"],
                ),
                PageBreak(),
            ]
        )

        story.extend(
            [
                Paragraph("2. Reproducible workflow and quality controls", styles["Heading1"]),
                Paragraph(
                    "The Nextflow workflow validates compressed-source MD5 values, performs read "
                    "QC/trimming, alignment, duplicate removal, MAPQ filtering, peak calling, "
                    "methylation extraction, DunedinPACE projection, SIRT6-locus occupancy "
                    "quantification, correlation analysis and manifest generation. The analysis "
                    "image pins tool/package versions and the DunedinPACE Git commit.",
                    styles["BodyText"],
                ),
                Spacer(1, 5 * mm),
                table(
                    [
                        ["Control", "Configured / measured value"],
                        ["Alignment MAPQ minimum", str(manifest.get("alignment", {}).get("mapq_threshold", "n/a"))],
                        ["ChIP peak q/FDR threshold", str(manifest.get("peak_calling", {}).get("fdr", "n/a"))],
                        ["SIRT6 CUT&RUN q/FDR threshold", str(manifest.get("cutrun", {}).get("fdr", "n/a"))],
                        ["Minimum CpG depth", str(pace_qc[0].get("min_cpg_depth", "n/a"))],
                    ],
                    [70 * mm, 98 * mm],
                ),
                PageBreak(),
            ]
        )

        story.extend(
            [
                Paragraph("3. DunedinPACE probe coverage QC", styles["Heading1"]),
                Paragraph(
                    "Each sample must retain sufficient observed background/model probes before "
                    "projection. Missing values are then handled by the official PACEProjector "
                    "implementation according to its documented cohort-level rules.",
                    styles["BodyText"],
                ),
                Spacer(1, 5 * mm),
                Image(str(qc_png), width=170 * mm, height=91 * mm),
                Spacer(1, 5 * mm),
                table(
                    [["Sample", "Matched", "Required", "Fraction"]]
                    + [
                        [
                            row["sample_id"],
                            row.get("matched_probes", ""),
                            row.get("required_probes", ""),
                            f"{finite_float(row['fraction'], 'probe fraction'):.3f}",
                        ]
                        for row in pace_qc
                    ],
                    [75 * mm, 30 * mm, 30 * mm, 30 * mm],
                ),
                PageBreak(),
            ]
        )

        story.extend(
            [
                Paragraph("4. ChIP-seq and SIRT6-locus analysis", styles["Heading1"]),
                Paragraph(
                    "H3K9ac and H3K56ac libraries are aligned to hg19, duplicate-filtered and "
                    "restricted to MAPQ >=30 before MACS3 peak calling against the condition-matched "
                    "INPUT. SIRT6-specific loci are independently derived from SIRT6 CUT&RUN "
                    "relative to IgG and then used as the occupancy measurement space. This avoids "
                    "choosing regions based on the desired DunedinPACE correlation.",
                    styles["BodyText"],
                ),
                Spacer(1, 7 * mm),
                Paragraph(
                    "Differential occupancy values shown in the following pages are computed by the "
                    "pipeline from aligned reads and called peaks; the report generator never alters "
                    "or rescales them to match an acceptance threshold.",
                    styles["BodyText"],
                ),
                PageBreak(),
            ]
        )

        story.extend(
            [
                Paragraph("5. Measured DunedinPACE scores", styles["Heading1"]),
                Image(str(scores_png), width=170 * mm, height=91 * mm),
                Spacer(1, 5 * mm),
                table(
                    [["Sample", "Condition", "Rep", "DunedinPACE"]]
                    + [
                        [
                            row["sample_id"],
                            row["condition"],
                            row["replicate"],
                            f"{finite_float(row['dunedinpace'], 'PACE'):.5f}",
                        ]
                        for row in scores
                    ],
                    [80 * mm, 35 * mm, 18 * mm, 35 * mm],
                ),
                PageBreak(),
            ]
        )

        for section, mark, image_path in (
            ("6", "H3K9ac", h3k9_png),
            ("7", "H3K56ac", h3k56_png),
        ):
            result = correlations[f"{mark}_vs_DunedinPACE"]
            r_value = finite_float(result["pearson_r"], f"{mark} r")
            p_value = finite_float(result["p_value"], f"{mark} p")
            n_value = int(result.get("n", correlations.get("n_paired", 0)))
            passes = r_value > 0.92 and p_value < 0.01
            story.extend(
                [
                    Paragraph(f"{section}. {mark} vs DunedinPACE", styles["Heading1"]),
                    Image(str(image_path), width=160 * mm, height=107 * mm),
                    Spacer(1, 5 * mm),
                    table(
                        [
                            ["Statistic", "Measured value"],
                            ["Paired observations", str(n_value)],
                            ["Pearson r", f"{r_value:.6f}"],
                            ["p-value", f"{p_value:.6g}"],
                            ["Requested r > 0.92 and p < 0.01", "PASS" if passes else "NOT MET"],
                        ],
                        [85 * mm, 83 * mm],
                    ),
                    PageBreak(),
                ]
            )

        story.extend(
            [
                Paragraph("8. Reproducibility, manifest and limitations", styles["Heading1"]),
                Paragraph(
                    "The machine-readable manifest is generated from the same measured correlation "
                    "and model-metadata files consumed by this report. A DOI and container checksum "
                    "are included only when real artifacts have been deposited/exported; the "
                    "pipeline intentionally does not invent placeholders that satisfy superficial "
                    "acceptance tests.",
                    styles["BodyText"],
                ),
                Spacer(1, 5 * mm),
                Paragraph(
                    "Interpretation limitations include tissue-domain mismatch for DunedinPACE, two "
                    "biological replicates per individual sirtuin condition in HRA003336, and the "
                    "fact that correlation does not establish causation. Results should therefore "
                    "be treated as an exploratory epigenomic association in an isogenic stem-cell "
                    "senescence model.",
                    styles["BodyText"],
                ),
                Spacer(1, 6 * mm),
                Paragraph("Manifest snapshot", styles["Heading2"]),
                Paragraph(
                    "<font name='Courier'>" + json.dumps(manifest, sort_keys=True)[:3500].replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;") + "</font>",
                    styles["Small"],
                ),
            ]
        )

        doc.build(story, onFirstPage=add_page_number, onLaterPages=add_page_number)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--occupancy", required=True, type=Path)
    parser.add_argument("--correlations", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--pace-qc", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    build_report(
        args.scores,
        args.occupancy,
        args.correlations,
        args.manifest,
        args.pace_qc,
        args.output,
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
