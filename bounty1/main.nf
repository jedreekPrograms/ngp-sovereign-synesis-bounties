nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: 'resources/samplesheet.csv'
params.controls_sheet = params.controls_sheet ?: 'resources/chip_inputs.csv'
params.cutrun_sheet = params.cutrun_sheet ?: 'resources/sirt6_cutrun.csv'
params.outdir = params.outdir ?: 'results'
params.reference_fasta_url = params.reference_fasta_url ?: 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz'
params.mapq = params.mapq ?: 30
params.peak_qvalue = params.peak_qvalue ?: 0.01
params.cutrun_qvalue = params.cutrun_qvalue ?: 0.05
params.min_cpg_depth = params.min_cpg_depth ?: 6
params.min_required_probe_fraction = params.min_required_probe_fraction ?: 0.80
params.sirt6_min_reciprocal_overlap = params.sirt6_min_reciprocal_overlap ?: 0.50
params.pace_probe_flank = params.pace_probe_flank ?: 500

process PREPARE_REFERENCE {
    tag 'hg19-reference'
    label 'big_mem'
    publishDir "${params.outdir}/provenance/reference", mode: 'copy'
    input:
    val reference_url
    output:
    path 'reference', emit: refdir
    script:
    """
    set -euo pipefail
    mkdir -p reference
    curl -fsSL --retry 5 --retry-delay 5 -o reference/hg19.fa.gz '${reference_url}'
    gzip -t reference/hg19.fa.gz
    gzip -d reference/hg19.fa.gz
    test -s reference/hg19.fa
    sha256sum reference/hg19.fa > reference/hg19.fa.sha256
    samtools faidx reference/hg19.fa
    bowtie2-build --threads ${task.cpus} reference/hg19.fa reference/hg19
    bwameth.py index reference/hg19.fa
    """
    stub:
    """
    mkdir -p reference
    printf '>chr1\nACGTACGTACGT\n' > reference/hg19.fa
    printf 'chr1\t12\t6\t12\t13\n' > reference/hg19.fa.fai
    echo 'stub' > reference/hg19.fa.sha256
    touch reference/hg19.1.bt2 reference/hg19.2.bt2 reference/hg19.3.bt2 reference/hg19.4.bt2
    touch reference/hg19.rev.1.bt2 reference/hg19.rev.2.bt2
    touch reference/hg19.fa.bwameth.c2t reference/hg19.fa.bwameth.c2t.amb
    """
}

process PREPARE_PACE_WINDOWS {
    tag 'DunedinPACE-probe-windows'
    publishDir "${params.outdir}/provenance", mode: 'copy'
    output:
    path 'dunedinpace_probe_windows.bed', emit: bed
    script:
    """
    set -euo pipefail
    Rscript ${projectDir}/bin/export_pace_regions.R \
      --output dunedinpace_probe_windows.bed \
      --flank ${params.pace_probe_flank}
    test -s dunedinpace_probe_windows.bed
    """
    stub:
    """
    printf 'chr1\t0\t1000\tcg00000029\n' > dunedinpace_probe_windows.bed
    """
}

process FASTP {
    tag "${meta.sample_id}"
    label 'stream_io'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', saveAs: { name ->
        (name.endsWith('.json') || name.endsWith('.html')) ? name : null
    }
    input:
    tuple val(meta), val(reads)
    output:
    tuple val(meta), path("${meta.sample_id}.trimmed*.fastq.gz"), emit: trimmed
    path "${meta.sample_id}.fastp.json", emit: json
    path "${meta.sample_id}.fastp.html", emit: html
    script:
    def r2 = reads.size() == 2 ? reads[1] : ''
    """
    set -euo pipefail
    bash ${projectDir}/bin/stream_fastp.sh \
      '${meta.sample_id}' \
      '${reads[0]}' \
      '${r2}' \
      '${meta.read1_md5}' \
      '${meta.read2_md5}' \
      '${task.cpus}'
    """
    stub:
    def outputs = meta.paired ? "touch ${meta.sample_id}.trimmed.R1.fastq.gz ${meta.sample_id}.trimmed.R2.fastq.gz" : "touch ${meta.sample_id}.trimmed.R1.fastq.gz"
    """
    ${outputs}
    echo '{}' > ${meta.sample_id}.fastp.json
    echo '<html></html>' > ${meta.sample_id}.fastp.html
    """
}

process ALIGN_CHIP {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/alignment/chip", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path refdir
    output:
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.bam"), path("${meta.sample_id}.mapq${params.mapq}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.flagstat.txt"), emit: flagstat
    script:
    def readArgs = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    """
    set -euo pipefail
    bowtie2 --very-sensitive -x ${refdir}/hg19 ${readArgs} -p ${task.cpus} 2> ${meta.sample_id}.bowtie2.log | \
      samtools view -b - | samtools sort -n -@ ${task.cpus} -o ${meta.sample_id}.name.bam
    samtools fixmate -m ${meta.sample_id}.name.bam ${meta.sample_id}.fixmate.bam
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}.position.bam ${meta.sample_id}.fixmate.bam
    samtools markdup -r -@ ${task.cpus} ${meta.sample_id}.position.bam ${meta.sample_id}.dedup.bam
    samtools view -b -q ${params.mapq} ${meta.sample_id}.dedup.bam | \
      samtools sort -@ ${task.cpus} -o ${meta.sample_id}.mapq${params.mapq}.bam
    samtools index ${meta.sample_id}.mapq${params.mapq}.bam
    samtools flagstat ${meta.sample_id}.mapq${params.mapq}.bam > ${meta.sample_id}.flagstat.txt
    rm -f ${meta.sample_id}.name.bam ${meta.sample_id}.fixmate.bam ${meta.sample_id}.position.bam ${meta.sample_id}.dedup.bam
    """
    stub:
    """
    touch ${meta.sample_id}.mapq${params.mapq}.bam
    touch ${meta.sample_id}.mapq${params.mapq}.bam.bai
    echo '0 + 0 mapped' > ${meta.sample_id}.flagstat.txt
    """
}

process CALL_PEAKS {
    tag "${meta.sample_id}:${meta.mark}"
    publishDir "${params.outdir}/peaks", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    path control_bams
    output:
    tuple val(meta), path("${meta.sample_id}_peaks.narrowPeak"), emit: peaks
    tuple val(meta), path("${meta.sample_id}_treat_pileup.bdg"), emit: pileup
    script:
    def macsFormat = meta.paired ? 'BAMPE' : 'BAM'
    def control = "${meta.condition}_INPUT.mapq${params.mapq}.bam"
    """
    set -euo pipefail
    test -f ${control}
    macs3 callpeak -t ${bam} -c ${control} -f ${macsFormat} -g hs \
      -n ${meta.sample_id} -q ${params.peak_qvalue} --keep-dup all -B
    """
    stub:
    """
    printf 'chr1\t100\t200\t${meta.sample_id}\t100\t.\t1\t1\t1\t50\n' > ${meta.sample_id}_peaks.narrowPeak
    printf 'chr1\t100\t200\t1\n' > ${meta.sample_id}_treat_pileup.bdg
    """
}

process ALIGN_CUTRUN {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/alignment/cutrun", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    path refdir
    output:
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.bam"), path("${meta.sample_id}.mapq${params.mapq}.bam.bai"), emit: bam
    script:
    def readArgs = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    """
    set -euo pipefail
    bowtie2 --very-sensitive --no-mixed --no-discordant --dovetail \
      -x ${refdir}/hg19 ${readArgs} -p ${task.cpus} 2> ${meta.sample_id}.bowtie2.log | \
      samtools view -b - | samtools sort -n -@ ${task.cpus} -o ${meta.sample_id}.name.bam
    samtools fixmate -m ${meta.sample_id}.name.bam ${meta.sample_id}.fixmate.bam
    samtools sort -@ ${task.cpus} -o ${meta.sample_id}.position.bam ${meta.sample_id}.fixmate.bam
    samtools markdup -r -@ ${task.cpus} ${meta.sample_id}.position.bam ${meta.sample_id}.dedup.bam
    samtools view -b -q ${params.mapq} ${meta.sample_id}.dedup.bam | \
      samtools sort -@ ${task.cpus} -o ${meta.sample_id}.mapq${params.mapq}.bam
    samtools index ${meta.sample_id}.mapq${params.mapq}.bam
    rm -f ${meta.sample_id}.name.bam ${meta.sample_id}.fixmate.bam ${meta.sample_id}.position.bam ${meta.sample_id}.dedup.bam
    """
    stub:
    """
    touch ${meta.sample_id}.mapq${params.mapq}.bam
    touch ${meta.sample_id}.mapq${params.mapq}.bam.bai
    """
}

process CALL_CUTRUN_PEAKS {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/cutrun/peaks", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    path control_bams
    output:
    tuple val(meta), path("${meta.sample_id}_peaks.narrowPeak"), emit: peaks
    script:
    def control = "${meta.condition}_rep${meta.replicate}_IgG_CUTRUN.mapq${params.mapq}.bam"
    """
    set -euo pipefail
    test -f ${control}
    macs3 callpeak -t ${bam} -c ${control} -f BAMPE -g hs \
      -n ${meta.sample_id} -q ${params.cutrun_qvalue} --keep-dup all --call-summits
    """
    stub:
    """
    printf 'chr1\t100\t200\t${meta.sample_id}\t100\t.\t1\t1\t1\t50\n' > ${meta.sample_id}_peaks.narrowPeak
    """
}

process BUILD_SIRT6_LOCI {
    tag 'SIRT6-specific-loci'
    publishDir "${params.outdir}/cutrun", mode: 'copy'
    input:
    path peak_files
    output:
    path 'sirt6_specific_loci.bed', emit: loci
    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/build_sirt6_loci.py \
      --peaks-glob '*_peaks.narrowPeak' \
      --min-reciprocal-overlap ${params.sirt6_min_reciprocal_overlap} \
      --output sirt6_specific_loci.bed
    """
    stub:
    """
    printf 'chr1\t120\t200\n' > sirt6_specific_loci.bed
    """
}

process ALIGN_WGBS_STREAM {
    tag "${meta.sample_id}"
    label 'stream_io'
    publishDir "${params.outdir}/alignment/wgbs", mode: 'copy'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', saveAs: { name ->
        (name.endsWith('.json') || name.endsWith('.html')) ? name : null
    }
    input:
    tuple val(meta), val(reads)
    path refdir
    path pace_bed
    output:
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.pace-targets.deduplicated.bam"), path("${meta.sample_id}.mapq${params.mapq}.pace-targets.deduplicated.bam.bai"), emit: bam
    path "${meta.sample_id}.fastp.json", emit: fastp_json
    path "${meta.sample_id}.fastp.html", emit: fastp_html
    path "${meta.sample_id}.pace-targets.flagstat.txt", emit: flagstat
    script:
    if (!meta.paired) {
        error "Streaming WGBS path requires paired-end input for ${meta.sample_id}"
    }
    """
    set -euo pipefail
    bash ${projectDir}/bin/stream_wgbs_bwameth.sh \
      '${meta.sample_id}' \
      '${reads[0]}' \
      '${reads[1]}' \
      '${meta.read1_md5}' \
      '${meta.read2_md5}' \
      '${refdir}/hg19.fa' \
      '${pace_bed}' \
      '${params.mapq}' \
      '${task.cpus}'
    """
    stub:
    """
    touch ${meta.sample_id}.mapq${params.mapq}.pace-targets.deduplicated.bam
    touch ${meta.sample_id}.mapq${params.mapq}.pace-targets.deduplicated.bam.bai
    echo '{}' > ${meta.sample_id}.fastp.json
    echo '<html></html>' > ${meta.sample_id}.fastp.html
    echo '0 + 0 mapped' > ${meta.sample_id}.pace-targets.flagstat.txt
    """
}

process EXTRACT_METHYLATION_STREAM {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/methylation", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    path refdir
    output:
    tuple val(meta), path("${meta.sample_id}.bismark.cov.gz"), emit: cov
    script:
    """
    set -euo pipefail
    micromamba run -n methyldackel MethylDackel extract \
      -q ${params.mapq} -p 5 --minDepth 1 -@ ${task.cpus} \
      ${refdir}/hg19.fa ${bam} -o ${meta.sample_id}.methyldackel
    test -s ${meta.sample_id}.methyldackel_CpG.bedGraph
    python3 ${projectDir}/bin/methyldackel_to_bismark.py \
      --input ${meta.sample_id}.methyldackel_CpG.bedGraph \
      --output ${meta.sample_id}.bismark.cov.gz
    gzip -t ${meta.sample_id}.bismark.cov.gz
    """
    stub:
    """
    printf 'chr1\t100\t100\t50\t3\t3\n' | gzip -c > ${meta.sample_id}.bismark.cov.gz
    """
}

process BUILD_PACE_MATRIX {
    tag 'DunedinPACE'
    label 'big_mem'
    publishDir "${params.outdir}/pace", mode: 'copy'
    input:
    path cov_files
    output:
    path 'dunedinpace_scores.csv', emit: scores
    path 'pace_probe_qc.csv', emit: qc
    path 'pace_model_metadata.csv', emit: model_metadata
    script:
    """
    set -euo pipefail
    Rscript ${projectDir}/bin/compute_dunedinpace.R \
      --cov-dir . \
      --scores dunedinpace_scores.csv \
      --qc pace_probe_qc.csv \
      --model-metadata pace_model_metadata.csv \
      --min-probe-fraction ${params.min_required_probe_fraction} \
      --min-depth ${params.min_cpg_depth}
    """
    stub:
    """
    cat > dunedinpace_scores.csv <<'EOF'
sample_id,condition,replicate,dunedinpace,model
WT_rep1_WGBS,WT,1,1.00,DunedinPACE
SIRT1_rep1_WGBS,SIRT1,1,1.10,DunedinPACE
SIRT2_rep1_WGBS,SIRT2,1,1.20,DunedinPACE
EOF
    cat > pace_probe_qc.csv <<'EOF'
sample_id,matched_probes,required_probes,fraction,min_cpg_depth
WT_rep1_WGBS,20000,20000,1.0,6
SIRT1_rep1_WGBS,20000,20000,1.0,6
SIRT2_rep1_WGBS,20000,20000,1.0,6
EOF
    cat > pace_model_metadata.csv <<'EOF'
model,intercept,required_background_probes,annotation,source_package
DunedinPACE,-1.949859,20000,IlluminaHumanMethylationEPICanno.ilm10b4.hg19,danbelsky/DunedinPACE@4b569983543e51d1022aecec9a25e694bb3a336a
EOF
    """
}

process CHIP_OCCUPANCY {
    tag 'SIRT6-locus-occupancy'
    publishDir "${params.outdir}/chip", mode: 'copy'
    input:
    path peak_files
    path bam_files
    path bai_files
    path loci_bed
    output:
    path 'histone_occupancy.csv', emit: occupancy
    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/compute_chip_occupancy.py \
      --peaks-glob '*_peaks.narrowPeak' \
      --bam-glob '*.bam' \
      --loci-bed ${loci_bed} \
      --output histone_occupancy.csv
    """
    stub:
    """
    cat > histone_occupancy.csv <<'EOF'
condition,replicate,mark,differential_occupancy,log2_cpm
WT,1,H3K9ac,0,1
WT,1,H3K56ac,0,1
SIRT1,1,H3K9ac,1,2
SIRT1,1,H3K56ac,1,2
SIRT2,1,H3K9ac,2,3
SIRT2,1,H3K56ac,2,3
EOF
    """
}

process CORRELATE {
    tag 'histone-vs-pace'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path scores
    path occupancy
    output:
    path 'correlations.json', emit: correlations
    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/correlate.py --pace ${scores} --occupancy ${occupancy} --output correlations.json
    """
    stub:
    """
    cat > correlations.json <<'EOF'
{"H3K9ac_vs_DunedinPACE":{"pearson_r":1.0,"p_value":0.0,"n":3},"H3K56ac_vs_DunedinPACE":{"pearson_r":1.0,"p_value":0.0,"n":3},"n_paired":3}
EOF
    """
}

process MANIFEST {
    tag 'manifest'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path correlations
    path model_metadata
    output:
    path 'manifest.json', emit: manifest
    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/build_manifest.py \
      --correlations ${correlations} \
      --model-metadata ${model_metadata} \
      --mapq ${params.mapq} \
      --peak-fdr ${params.peak_qvalue} \
      --cutrun-fdr ${params.cutrun_qvalue} \
      --sirt6-min-reciprocal-overlap ${params.sirt6_min_reciprocal_overlap} \
      --output manifest.json
    """
    stub:
    """
    cat > manifest.json <<'EOF'
{"alignment":{"mapq_threshold":30},"peak_calling":{"fdr":0.01},"cutrun":{"fdr":0.05}}
EOF
    """
}

process REPORT {
    tag 'supplementary-report'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path scores
    path occupancy
    path correlations
    path manifest
    path pace_qc
    output:
    path 'supplementary_report.pdf', emit: pdf
    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/generate_report.py \
      --scores ${scores} \
      --occupancy ${occupancy} \
      --correlations ${correlations} \
      --manifest ${manifest} \
      --pace-qc ${pace_qc} \
      --output supplementary_report.pdf
    test -s supplementary_report.pdf
    """
    stub:
    """
    printf '%s' '%PDF-1.4 stub' > supplementary_report.pdf
    """
}

workflow {
    primary = Channel.fromPath(params.samplesheet, checkIfExists: true).splitCsv(header: true).map { row ->
        def meta = [
            sample_id: row.sample_id,
            condition: row.condition,
            replicate: row.replicate as Integer,
            assay: row.assay,
            mark: row.mark ?: 'NA',
            paired: row.fastq_2 ? true : false,
            read1_md5: row.read1_md5 ?: '',
            read2_md5: row.read2_md5 ?: ''
        ]
        def reads = row.fastq_2 ? [row.fastq_1, row.fastq_2] : [row.fastq_1]
        tuple(meta, reads)
    }

    controls = Channel.fromPath(params.controls_sheet, checkIfExists: true).splitCsv(header: true).map { row ->
        def meta = [
            sample_id: row.sample_id,
            condition: row.condition,
            replicate: row.replicate as Integer,
            assay: row.assay,
            mark: 'NA',
            paired: row.fastq_2 ? true : false,
            read1_md5: row.read1_md5 ?: '',
            read2_md5: row.read2_md5 ?: ''
        ]
        def reads = row.fastq_2 ? [row.fastq_1, row.fastq_2] : [row.fastq_1]
        tuple(meta, reads)
    }

    cutrun = Channel.fromPath(params.cutrun_sheet, checkIfExists: true).splitCsv(header: true).map { row ->
        def meta = [
            sample_id: row.sample_id,
            condition: row.condition,
            replicate: row.replicate as Integer,
            assay: row.assay,
            target: row.target,
            mark: row.target,
            control: row.control?.toString()?.toLowerCase() == 'true',
            paired: row.fastq_2 ? true : false,
            read1_md5: row.read1_md5 ?: '',
            read2_md5: row.read2_md5 ?: ''
        ]
        def reads = row.fastq_2 ? [row.fastq_1, row.fastq_2] : [row.fastq_1]
        tuple(meta, reads)
    }

    ref_queue = PREPARE_REFERENCE(Channel.value(params.reference_fasta_url)).refdir
    ref_value = ref_queue.collect().map { items -> items[0] }
    pace_queue = PREPARE_PACE_WINDOWS().bed
    pace_value = pace_queue.collect().map { items -> items[0] }

    wgbs_raw = primary.filter { meta, reads -> meta.assay == 'WGBS' }
    chip_primary = primary.filter { meta, reads -> meta.assay == 'CHIP' }
    non_wgbs = chip_primary.mix(controls).mix(cutrun)
    trimmed = FASTP(non_wgbs).trimmed

    chip_and_input_trim = trimmed.filter { meta, reads -> meta.assay == 'CHIP' || meta.assay == 'INPUT' }
    cutrun_trim = trimmed.filter { meta, reads -> meta.assay == 'CUTRUN' }

    chip_all_bam = ALIGN_CHIP(chip_and_input_trim, ref_value).bam
    chip_bam = chip_all_bam.filter { meta, bam, bai -> meta.assay == 'CHIP' }
    input_bam = chip_all_bam.filter { meta, bam, bai -> meta.assay == 'INPUT' }
    control_bams = input_bam.map { meta, bam, bai -> bam }.collect()
    peaks = CALL_PEAKS(chip_bam, control_bams).peaks

    cutrun_all_bam = ALIGN_CUTRUN(cutrun_trim, ref_value).bam
    cutrun_target_bam = cutrun_all_bam.filter { meta, bam, bai -> meta.target == 'SIRT6' }
    cutrun_igg_bam = cutrun_all_bam.filter { meta, bam, bai -> meta.target == 'IgG' }
    cutrun_controls = cutrun_igg_bam.map { meta, bam, bai -> bam }.collect()
    cutrun_peaks = CALL_CUTRUN_PEAKS(cutrun_target_bam, cutrun_controls).peaks
    sirt6_loci = BUILD_SIRT6_LOCI(cutrun_peaks.map { meta, peak -> peak }.collect()).loci

    wgbs_target_bam = ALIGN_WGBS_STREAM(wgbs_raw, ref_value, pace_value).bam
    cov = EXTRACT_METHYLATION_STREAM(wgbs_target_bam, ref_value).cov
    pace_result = BUILD_PACE_MATRIX(cov.map { meta, coverage -> coverage }.collect())

    occupancy = CHIP_OCCUPANCY(
        peaks.map { meta, peak -> peak }.collect(),
        chip_bam.map { meta, bam, bai -> bam }.collect(),
        chip_bam.map { meta, bam, bai -> bai }.collect(),
        sirt6_loci
    ).occupancy

    corr = CORRELATE(pace_result.scores, occupancy).correlations
    manifest = MANIFEST(corr, pace_result.model_metadata).manifest
    REPORT(pace_result.scores, occupancy, corr, manifest, pace_result.qc)
}
