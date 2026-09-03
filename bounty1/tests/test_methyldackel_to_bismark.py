import gzip
import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "bin" / "methyldackel_to_bismark.py"
SPEC = importlib.util.spec_from_file_location("methyldackel_to_bismark", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_convert_line_changes_zero_based_to_one_based():
    row = MODULE.convert_line("chr22\t99\t100\t75\t3\t1\n")
    assert row == "chr22\t100\t100\t75\t3\t1\n"


def test_convert_line_skips_metadata():
    assert MODULE.convert_line("track type=bedGraph\n") is None
    assert MODULE.convert_line("# comment\n") is None
    assert MODULE.convert_line("\n") is None


def test_convert_file_writes_gzip_coverage(tmp_path):
    source = tmp_path / "sample_CpG.bedGraph"
    source.write_text(
        "track type=bedGraph\nchr1\t9\t10\t50\t2\t2\nchr2\t19\t20\t100\t6\t0\n",
        encoding="utf-8",
    )
    output = tmp_path / "sample.bismark.cov.gz"
    rows = MODULE.convert_file(source, output)
    assert rows == 2
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        assert handle.read() == "chr1\t10\t10\t50\t2\t2\nchr2\t20\t20\t100\t6\t0\n"
