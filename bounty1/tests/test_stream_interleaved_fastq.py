import gzip
import hashlib
import importlib.util
import io
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'bin' / 'stream_interleaved_fastq.py'
spec = importlib.util.spec_from_file_location('stream_interleaved_fastq', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def write_fastq(path, mate, count=3):
    with gzip.open(path, 'wb') as handle:
        for index in range(count):
            sequence = b'ACGTACGT' if mate == 1 else b'TGCATGCA'
            handle.write(f'@read{index}/{mate}\n'.encode())
            handle.write(sequence + b'\n+\n' + (b'I' * len(sequence)) + b'\n')


def md5(path):
    digest = hashlib.md5(usedforsecurity=False)
    digest.update(Path(path).read_bytes())
    return digest.hexdigest()


def test_paired_stream_is_interleaved_and_validates_compressed_md5(tmp_path):
    r1 = tmp_path / 'R1.fq.gz'
    r2 = tmp_path / 'R2.fq.gz'
    write_fastq(r1, 1)
    write_fastq(r2, 2)
    out = io.BytesIO()

    count = module.stream_fastq(
        str(r1), str(r2), out,
        r1_md5=md5(r1), r2_md5=md5(r2), timeout=5,
    )

    assert count == 3
    lines = out.getvalue().splitlines()
    assert lines[0] == b'@read0/1'
    assert lines[4] == b'@read0/2'
    assert lines[8] == b'@read1/1'
    assert lines[12] == b'@read1/2'


def test_stream_rejects_archive_md5_mismatch(tmp_path):
    r1 = tmp_path / 'R1.fq.gz'
    write_fastq(r1, 1, count=1)

    try:
        module.stream_fastq(str(r1), '', io.BytesIO(), r1_md5='0' * 32, timeout=5)
        assert False, 'expected an MD5 mismatch'
    except ValueError as exc:
        assert 'MD5 mismatch' in str(exc)


def test_stream_rejects_mismatched_pair_names(tmp_path):
    r1 = tmp_path / 'R1.fq.gz'
    r2 = tmp_path / 'R2.fq.gz'
    write_fastq(r1, 1, count=1)
    with gzip.open(r2, 'wb') as handle:
        handle.write(b'@other/2\nTGCATGCA\n+\nIIIIIIII\n')

    try:
        module.stream_fastq(str(r1), str(r2), io.BytesIO(), timeout=5)
        assert False, 'expected pair-name mismatch'
    except ValueError as exc:
        assert 'pair-name mismatch' in str(exc)
