import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import correlate


def test_pearson_requires_three():
    try:
        correlate._pearson([1, 2], [1, 2])
        assert False
    except ValueError:
        pass


def test_compute_pairs_by_condition_and_replicate():
    pace = {('WT',1):0.9, ('SIRT1',1):1.0, ('SIRT2',1):1.1, ('SIRT3',1):1.2}
    occ = {
        ('WT',1): {'H3K9ac':0.0,'H3K56ac':0.0},
        ('SIRT1',1): {'H3K9ac':1.0,'H3K56ac':1.1},
        ('SIRT2',1): {'H3K9ac':2.0,'H3K56ac':2.1},
        ('SIRT3',1): {'H3K9ac':3.0,'H3K56ac':3.1},
    }
    out = correlate.compute(pace, occ)
    assert out['n_paired'] == 4
    assert out['H3K9ac_vs_DunedinPACE']['pearson_r'] > 0.99


def test_manifest_marks_data_as_non_synthetic(tmp_path):
    c = tmp_path/'c.json'
    c.write_text(json.dumps({
        'n_paired': 4,
        'H3K9ac_vs_DunedinPACE': {'pearson_r': 0.95, 'p_value': 0.005},
        'H3K56ac_vs_DunedinPACE': {'pearson_r': 0.96, 'p_value': 0.004}
    }))
    o = tmp_path/'manifest.json'
    subprocess.check_call([
        sys.executable, str(ROOT/'bin/build_manifest.py'),
        '--correlations', str(c), '--mapq','30','--peak-fdr','0.01','--output',str(o)
    ])
    m = json.loads(o.read_text())
    assert m['alignment']['mapq_threshold'] == 30
    assert m['peak_calling']['fdr'] < 0.05
    assert m['provenance']['synthetic_values_used'] is False
    assert m['data_deposit_doi'] == ''


def test_r_pipeline_uses_upstream_api_and_background_probes():
    text = (ROOT/'bin/compute_dunedinpace.R').read_text()
    assert 'PACEProjector(' in text
    assert 'getRequiredProbes(backgroundList=TRUE)' in text
    assert 'DunedinPACE.getRequiredProbes' not in text
