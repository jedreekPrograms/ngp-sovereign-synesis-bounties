#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from pathlib import Path


def load_model_metadata(path):
    with Path(path).open(newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 1:
        raise ValueError(f'expected exactly one DunedinPACE model metadata row, found {len(rows)}')
    row = rows[0]
    intercept = float(row['intercept'])
    return {
        'model': row['model'],
        'intercept': intercept,
        'required_background_probes': int(row['required_background_probes']),
        'annotation': row['annotation'],
        'source_package': row['source_package'],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--correlations', required=True)
    ap.add_argument('--model-metadata', required=True)
    ap.add_argument('--mapq', required=True, type=int)
    ap.add_argument('--peak-fdr', required=True, type=float)
    ap.add_argument('--cutrun-fdr', default=0.05, type=float)
    ap.add_argument('--sirt6-min-reciprocal-overlap', default=0.50, type=float)
    ap.add_argument('--doi', default='')
    ap.add_argument('--docker-image', default='')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()

    corr = json.loads(Path(a.correlations).read_text(encoding='utf-8'))
    pace_model = load_model_metadata(a.model_metadata)

    manifest = {
        'pipeline_version': '1.0.0',
        'data_source': {
            'primary_project': 'PRJCA012536',
            'primary_study': 'HRA003336',
            'sirt6_cutrun_study': 'HRA005392',
            'access': 'open',
            'description': (
                'paired H3K9ac/H3K56ac ChIP-seq and WGBS in WT and '
                'SIRT1-7-deficient human mesenchymal progenitor cells, with '
                'independent SIRT6 CUT&RUN used to define occupancy loci'
            ),
        },
        'locus_definition': {
            'assay': 'SIRT6 CUT&RUN',
            'study': 'HRA005392',
            'strategy': (
                'reproducible WT SIRT6 peaks across two replicates with '
                'reproducible SIRT6-KO antibody peaks removed'
            ),
            'minimum_reciprocal_overlap': a.sirt6_min_reciprocal_overlap,
            'independent_of_histone_marks': True,
        },
        'dunedinpace': {
            'intercept': pace_model['intercept'],
            'model': pace_model['model'],
            'implementation': pace_model['source_package'],
            'annotation': pace_model['annotation'],
            'required_background_probes': pace_model['required_background_probes'],
            'input': 'WGBS-derived beta values at required CpG probes',
        },
        'correlations': {
            'H3K9ac_vs_DunedinPACE': corr['H3K9ac_vs_DunedinPACE'],
            'H3K56ac_vs_DunedinPACE': corr['H3K56ac_vs_DunedinPACE'],
        },
        'n_paired_samples': corr['n_paired'],
        'peak_calling': {
            'caller': 'macs3',
            'fdr': a.peak_fdr,
            'qvalue_threshold': a.peak_fdr,
            'sirt6_cutrun_qvalue_threshold': a.cutrun_fdr,
        },
        'alignment': {
            'aligner': 'bowtie2/bismark',
            'reference_genome': 'hg19',
            'mapq_threshold': a.mapq,
            'duplicate_policy': 'PCR duplicates removed before MAPQ filtering for ChIP-seq and CUT&RUN; Bismark deduplication for WGBS',
        },
        'data_deposit_doi': a.doi,
        'provenance': {
            'correlations_computed': True,
            'dunedinpace_intercept_read_from_upstream_model': True,
            'sirt6_loci_derived_independently': True,
            'synthetic_values_used': False,
        },
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
