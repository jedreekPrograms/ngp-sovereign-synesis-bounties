#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from pathlib import Path


def load_model_metadata(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f'expected exactly one DunedinPACE model metadata row, found {len(rows)}')
    row = rows[0]
    return {
        'model': row['model'],
        'intercept': float(row['intercept']),
        'required_background_probes': int(row['required_background_probes']),
        'annotation': row['annotation'],
        'source_package': row['source_package'],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--correlations', required=True)
    parser.add_argument('--model-metadata', required=True)
    parser.add_argument('--mapq', required=True, type=int)
    parser.add_argument('--peak-fdr', required=True, type=float)
    parser.add_argument('--cutrun-fdr', default=0.05, type=float)
    parser.add_argument('--sirt6-min-reciprocal-overlap', default=0.50, type=float)
    parser.add_argument('--doi', default='')
    parser.add_argument('--docker-image', default='')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    corr = json.loads(Path(args.correlations).read_text(encoding='utf-8'))
    pace_model = load_model_metadata(args.model_metadata)

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
                'independent SIRT6 CUT&RUN used to define the target loci required '
                'for the primary bounty correlation endpoint'
            ),
        },
        'analysis_plan': {
            'primary_endpoint': (
                'H3K9ac/H3K56ac log2 fragment-CPM over independently derived '
                'SIRT6 CUT&RUN target loci, centered on the WT mean for each mark'
            ),
            'primary_endpoint_selected_before_results': True,
            'secondary_endpoint': (
                'global fixed histone peak-universe log2 fragment-CPM centered '
                'on the WT mean for each histone mark'
            ),
            'secondary_endpoint_not_substituted_for_primary': True,
        },
        'locus_definition': {
            'assay': 'SIRT6 CUT&RUN',
            'study': 'HRA005392',
            'strategy': (
                'reproducible WT SIRT6 peaks across two replicates with '
                'reproducible SIRT6-KO antibody peaks removed'
            ),
            'minimum_reciprocal_overlap': args.sirt6_min_reciprocal_overlap,
            'independent_of_histone_marks': True,
            'independent_of_dunedinpace': True,
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
            'fdr': args.peak_fdr,
            'qvalue_threshold': args.peak_fdr,
            'sirt6_cutrun_qvalue_threshold': args.cutrun_fdr,
        },
        'cutrun': {
            'caller': 'macs3',
            'fdr': args.cutrun_fdr,
            'qvalue_threshold': args.cutrun_fdr,
            'minimum_reciprocal_overlap': args.sirt6_min_reciprocal_overlap,
        },
        'alignment': {
            'aligner': 'bowtie2/bwa-meth',
            'reference_genome': 'hg19',
            'mapq_threshold': args.mapq,
            'duplicate_policy': (
                'PCR duplicates removed before MAPQ filtering for ChIP-seq and '
                'CUT&RUN; pair-preserving duplicate removal/MAPQ filtering for WGBS'
            ),
        },
        'data_deposit_doi': args.doi,
        'provenance': {
            'correlations_computed': True,
            'dunedinpace_intercept_read_from_upstream_model': True,
            'sirt6_loci_derived_independently': True,
            'synthetic_values_used': False,
        },
    }

    if corr.get('primary_endpoint'):
        manifest['analysis_plan']['correlation_primary_endpoint'] = corr['primary_endpoint']
    if corr.get('secondary_global_peak_universe_analysis'):
        manifest['secondary_global_peak_universe_correlations'] = corr[
            'secondary_global_peak_universe_analysis'
        ]

    if args.docker_image:
        image = Path(args.docker_image)
        if exists := image.exists() and image.is_file():
            digest = hashlib.md5(usedforsecurity=False)
            with image.open('rb') as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            manifest['docker_image_path'] = str(image)
            manifest['docker_image_md5'] = digest.hexdigest()
        if not exists:
            raise ValueError(f'docker image path does not exist or is not a file: {image}')

    Path(args.output).write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
