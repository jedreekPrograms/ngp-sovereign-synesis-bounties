import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import correlate
import compute_chip_occupancy as occupancy


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


def test_fragment_representative_counts_one_mate_per_pair():
    class Read:
        def __init__(self, **kw):
            self.is_unmapped = kw.get('is_unmapped', False)
            self.is_secondary = kw.get('is_secondary', False)
            self.is_supplementary = kw.get('is_supplementary', False)
            self.is_paired = kw.get('is_paired', True)
            self.is_proper_pair = kw.get('is_proper_pair', True)
            self.is_read1 = kw.get('is_read1', True)

    assert occupancy.is_fragment_representative(Read(is_read1=True)) is True
    assert occupancy.is_fragment_representative(Read(is_read1=False)) is False
    assert occupancy.is_fragment_representative(Read(is_proper_pair=False)) is False
    assert occupancy.is_fragment_representative(Read(is_secondary=True)) is False
    assert occupancy.is_fragment_representative(Read(is_supplementary=True)) is False
    assert occupancy.is_fragment_representative(Read(is_paired=False)) is True


def test_manifest_uses_model_metadata_and_marks_data_as_non_synthetic(tmp_path):
    c = tmp_path/'c.json'
    c.write_text(json.dumps({
        'n_paired': 4,
        'H3K9ac_vs_DunedinPACE': {'pearson_r': 0.95, 'p_value': 0.005},
        'H3K56ac_vs_DunedinPACE': {'pearson_r': 0.96, 'p_value': 0.004}
    }))
    model = tmp_path/'model.csv'
    with model.open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=['model','intercept','required_background_probes','annotation','source_package'],
        )
        w.writeheader()
        w.writerow({
            'model': 'fixture-model',
            'intercept': '-1.2345',
            'required_background_probes': '20000',
            'annotation': 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19',
            'source_package': 'danbelsky/DunedinPACE',
        })
    o = tmp_path/'manifest.json'
    subprocess.check_call([
        sys.executable, str(ROOT/'bin/build_manifest.py'),
        '--correlations', str(c), '--model-metadata', str(model),
        '--mapq','30','--peak-fdr','0.01','--output',str(o)
    ])
    m = json.loads(o.read_text())
    assert m['alignment']['mapq_threshold'] == 30
    assert m['alignment']['reference_genome'] == 'hg19'
    assert m['peak_calling']['fdr'] < 0.05
    assert m['dunedinpace']['intercept'] == -1.2345
    assert m['dunedinpace']['model'] == 'fixture-model'
    assert m['data_source']['primary_study'] == 'HRA003336'
    assert m['data_source']['sirt6_cutrun_study'] == 'HRA005392'
    assert m['locus_definition']['independent_of_histone_marks'] is True
    assert m['locus_definition']['minimum_reciprocal_overlap'] == 0.50
    assert m['provenance']['dunedinpace_intercept_read_from_upstream_model'] is True
    assert m['provenance']['sirt6_loci_derived_independently'] is True
    assert m['provenance']['synthetic_values_used'] is False
    assert m['data_deposit_doi'] == ''


def test_r_pipeline_uses_upstream_api_background_probes_and_model_intercept():
    text = (ROOT/'bin/compute_dunedinpace.R').read_text()
    assert 'PACEProjector(' in text
    assert 'getRequiredProbes(backgroundList=TRUE)' in text
    assert 'mPACE_Models$model_intercept[[model_name]]' in text
    assert 'DunedinPACE.getRequiredProbes' not in text


def test_paired_wgbs_stays_name_sorted_for_methylation_extraction():
    text = (ROOT/'main.nf').read_text()
    wgbs = text.split('process ALIGN_WGBS {', 1)[1].split('process EXTRACT_METHYLATION {', 1)[0]
    assert 'samtools sort -n' in wgbs
    assert 'deduplicated.name.bam' in wgbs
    assert 'samtools index' not in wgbs
    extractor = text.split('process EXTRACT_METHYLATION {', 1)[1].split('process BUILD_PACE_MATRIX {', 1)[0]
    assert '--paired-end' in extractor
