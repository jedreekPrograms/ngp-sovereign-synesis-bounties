import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ACCESSIONS = ROOT / "resources" / "hra003336_accessions.csv"


def test_hra003336_accession_manifest_is_complete_and_unique():
    with ACCESSIONS.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    assert len(rows) == 48
    assert len({row["sample_id"] for row in rows}) == 48
    assert len({row["biosample_accession"] for row in rows}) == 48
    assert all(row["study_accession"] == "HRA003336" for row in rows)
    assert all(row["bioproject_accession"] == "PRJCA012536" for row in rows)

    paired = defaultdict(set)
    for row in rows:
        key = (row["condition"], int(row["replicate"]))
        if row["assay"] == "WGBS":
            paired[key].add("WGBS")
        else:
            paired[key].add(row["mark"])

    expected = {"H3K9ac", "H3K56ac", "WGBS"}
    assert len(paired) == 16
    assert all(assays == expected for assays in paired.values())
