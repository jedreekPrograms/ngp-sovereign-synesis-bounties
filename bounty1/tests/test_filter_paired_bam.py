import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin' / 'filter_paired_bam.py'
spec = importlib.util.spec_from_file_location('filter_paired_bam', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class Read:
    def __init__(
        self,
        *,
        read1=False,
        read2=False,
        mapq=60,
        unmapped=False,
        secondary=False,
        supplementary=False,
    ):
        self.is_read1 = read1
        self.is_read2 = read2
        self.mapping_quality = mapq
        self.is_unmapped = unmapped
        self.is_secondary = secondary
        self.is_supplementary = supplementary


def test_pair_is_kept_only_when_both_mates_pass_mapq():
    assert module.keep_paired_group(
        [Read(read1=True, mapq=30), Read(read2=True, mapq=42)], 30
    ) is True
    assert module.keep_paired_group(
        [Read(read1=True, mapq=30), Read(read2=True, mapq=29)], 30
    ) is False


def test_pair_rejects_orphans_unmapped_and_non_primary_alignments():
    assert module.keep_paired_group([Read(read1=True)], 30) is False
    assert module.keep_paired_group(
        [Read(read1=True), Read(read2=True, unmapped=True)], 30
    ) is False
    assert module.keep_paired_group(
        [Read(read1=True), Read(read2=True), Read(secondary=True)], 30
    ) is True
    assert module.keep_paired_group(
        [Read(read1=True), Read(read2=True), Read(supplementary=True)], 30
    ) is True


def test_pair_requires_one_r1_and_one_r2():
    assert module.keep_paired_group(
        [Read(read1=True), Read(read1=True)], 30
    ) is False
    assert module.keep_paired_group(
        [Read(read2=True), Read(read2=True)], 30
    ) is False
