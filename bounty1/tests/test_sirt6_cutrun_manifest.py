import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "resources" / "sirt6_cutrun.csv"


def _rows():
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_sirt6_cutrun_manifest_has_expected_design():
    rows = _rows()
    assert len(rows) == 8
    assert {row["condition"] for row in rows} == {"WT", "SIRT6_KO"}
    assert {int(row["replicate"]) for row in rows} == {1, 2}
    assert {row["target"] for row in rows} == {"SIRT6", "IgG"}

    keys = {(row["condition"], int(row["replicate"]), row["target"]) for row in rows}
    expected = {
        (condition, replicate, target)
        for condition in ("WT", "SIRT6_KO")
        for replicate in (1, 2)
        for target in ("SIRT6", "IgG")
    }
    assert keys == expected


def test_sirt6_cutrun_accessions_and_md5_are_unique_and_well_formed():
    rows = _rows()
    for field, prefix in (
        ("biosample_accession", "HRS"),
        ("experiment_accession", "HRX"),
        ("run_accession", "HRR"),
    ):
        values = [row[field] for row in rows]
        assert len(values) == len(set(values))
        assert all(value.startswith(prefix) for value in values)

    md5_values = [row[field] for row in rows for field in ("read1_md5", "read2_md5")]
    assert all(len(value) == 32 and all(c in "0123456789abcdef" for c in value) for value in md5_values)


def test_sirt6_cutrun_r2_urls_match_r2_filenames():
    # HRA005392's exported metadata duplicates the R1 URL in the R2 URL column.
    # The repository manifest intentionally reconstructs R2 URLs from the run
    # accession plus the authoritative R2 filename and keeps the archive MD5.
    for row in _rows():
        run = row["run_accession"]
        assert row["fastq_1"].endswith(f"/{run}_f1.fq.gz")
        assert row["fastq_2"].endswith(f"/{run}_r2.fq.gz")
        assert row["source_study"] == "HRA005392"
        assert row["assay"] == "CUTRUN"
        assert (row["target"] == "IgG") == (row["control"].lower() == "true")
