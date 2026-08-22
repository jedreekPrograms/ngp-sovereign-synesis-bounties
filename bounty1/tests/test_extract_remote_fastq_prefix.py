import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin' / 'extract_remote_fastq_prefix.py'
spec = importlib.util.spec_from_file_location('extract_remote_fastq_prefix', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write_fastq(path, mate, count=5):
    with gzip.open(path, 'wb') as handle:
        for i in range(count):
            name = f'@read{i}/{mate}\n'.encode()
            seq = (b'ACGT' if mate == 1 else b'TGCA') + b'ACGT\n'
            qual = b'I' * (len(seq.rstrip())) + b'\n'
            handle.write(name)
            handle.write(seq)
            handle.write(b'+\n')
            handle.write(qual)


def _count_records(path):
    with gzip.open(path, 'rt') as handle:
        return sum(1 for _ in handle) // 4


def test_extract_local_file_urls_preserves_pairs_and_limit(tmp_path):
    src1 = tmp_path / 'source_R1.fastq.gz'
    src2 = tmp_path / 'source_R2.fastq.gz'
    out1 = tmp_path / 'out_R1.fastq.gz'
    out2 = tmp_path / 'out_R2.fastq.gz'
    _write_fastq(src1, 1)
    _write_fastq(src2, 2)

    module.extract(src1.as_uri(), src2.as_uri(), out1, out2, reads=3, timeout=10)

    assert _count_records(out1) == 3
    assert _count_records(out2) == 3
    with gzip.open(out1, 'rt') as h1, gzip.open(out2, 'rt') as h2:
        names1 = [h1.readline().strip() for _ in range(3) for __ in [0] if not [h1.readline(), h1.readline(), h1.readline()]]
        # The compact comprehension above is intentionally not used for pair
        # comparison; reopen and read records explicitly for clarity below.
    with gzip.open(out1, 'rt') as h1, gzip.open(out2, 'rt') as h2:
        for i in range(3):
            r1 = [h1.readline() for _ in range(4)]
            r2 = [h2.readline() for _ in range(4)]
            assert module.normalise_read_name(r1[0].encode()) == module.normalise_read_name(r2[0].encode())


def test_extract_rejects_mismatched_pairs(tmp_path):
    src1 = tmp_path / 'source_R1.fastq.gz'
    src2 = tmp_path / 'source_R2.fastq.gz'
    out1 = tmp_path / 'out_R1.fastq.gz'
    out2 = tmp_path / 'out_R2.fastq.gz'
    _write_fastq(src1, 1, count=1)
    with gzip.open(src2, 'wb') as handle:
        handle.write(b'@different/2\nTGCAACGT\n+\nIIIIIIII\n')

    try:
        module.extract(src1.as_uri(), src2.as_uri(), out1, out2, reads=1, timeout=10)
        assert False, 'expected pair mismatch'
    except ValueError as exc:
        assert 'pair-name mismatch' in str(exc)
