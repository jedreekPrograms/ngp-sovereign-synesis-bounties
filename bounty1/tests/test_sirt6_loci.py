import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import build_sirt6_loci


def _write_peak(path, rows):
    with path.open('w', encoding='utf-8') as handle:
        for chrom, start, end in rows:
            handle.write(f'{chrom}\t{start}\t{end}\tpeak\t100\t.\t1\t1\t1\t1\n')


def test_build_loci_requires_reproducible_wt_and_removes_reproducible_ko(tmp_path):
    files = []
    definitions = {
        'WT_rep1_SIRT6_CUTRUN_peaks.narrowPeak': [
            ('chr1', 100, 200),
            ('chr1', 500, 600),
            ('chr2', 1000, 1100),
        ],
        'WT_rep2_SIRT6_CUTRUN_peaks.narrowPeak': [
            ('chr1', 120, 220),
            ('chr1', 510, 590),
            ('chr2', 1020, 1120),
        ],
        'SIRT6_KO_rep1_SIRT6_CUTRUN_peaks.narrowPeak': [
            ('chr1', 505, 585),
        ],
        'SIRT6_KO_rep2_SIRT6_CUTRUN_peaks.narrowPeak': [
            ('chr1', 515, 575),
        ],
    }
    for name, rows in definitions.items():
        path = tmp_path / name
        _write_peak(path, rows)
        files.append(str(path))

    loci = build_sirt6_loci.build_loci(files, min_reciprocal=0.50)

    assert ('chr1', 120, 200) in loci
    assert ('chr2', 1020, 1100) in loci
    assert all(not (chrom == 'chr1' and start < 590 and end > 510) for chrom, start, end in loci)


def test_reciprocal_overlap_rejects_tiny_edge_overlap():
    a = ('chr1', 100, 200)
    b = ('chr1', 190, 290)
    assert build_sirt6_loci.reciprocal_overlap(a, b) == 0.1
    assert build_sirt6_loci.reproducible_intersections([a], [b], 0.50) == []
