#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

DUNEDINPACE_INTERCEPT = 51.024577


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--correlations', required=True)
    ap.add_argument('--mapq', required=True, type=int)
    ap.add_argument('--peak-fdr', required=True, type=float)
    ap.add_argument('--doi', default='')
    ap.add_argument('--docker-image', default='')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    corr = json.loads(Path(a.correlations).read_text(encoding='utf-8'))
    manifest = {
        'pipeline_version': '1.0.0',
        'data_source': {
            'project': 'PRJCA012536',
            'study': 'HRA003336',
            'access': 'open',
            'description': 'paired H3K9ac/H3K56ac ChIP-seq and WGBS in WT and SIRT1-7-deficient human mesenchymal progenitor cells'
        },
        'dunedinpace': {
            'intercept': DUNEDINPACE_INTERCEPT,
            'implementation': 'danbelsky/DunedinPACE',
            'input': 'WGBS-derived beta values at required CpG probes'
        },
        'correlations': {
            'H3K9ac_vs_DunedinPACE': corr['H3K9ac_vs_DunedinPACE'],
            'H3K56ac_vs_DunedinPACE': corr['H3K56ac_vs_DunedinPACE']
        },
        'n_paired_samples': corr['n_paired'],
        'peak_calling': {'caller': 'macs3', 'fdr': a.peak_fdr, 'qvalue_threshold': a.peak_fdr},
        'alignment': {'aligner': 'bowtie2', 'mapq_threshold': a.mapq},
        'data_deposit_doi': a.doi,
        'provenance': {
            'correlations_computed': True,
            'synthetic_values_used': False
        }
    }
    if a.docker_image:
        p = Path(a.docker_image)
        if p.exists() and p.is_file():
            h = hashlib.md5()
            with p.open('rb') as fh:
                for chunk in iter(lambda: fh.read(1024 * 1024), b''):
                    h.update(chunk)
            manifest['docker_image_path'] = str(p)
            manifest['docker_image_md5'] = h.hexdigest()
    Path(a.output).write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
