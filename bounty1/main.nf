nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: 'resources/samplesheet.csv'
params.controls_sheet = params.controls_sheet ?: 'resources/chip_inputs.csv'
params.cutrun_sheet = params.cutrun_sheet ?: 'resources/sirt6_cutrun.csv'
params.outdir = params.outdir ?: 'results'
params.mapq = params.mapq ?: 30
params.peak_qvalue = params.peak_qvalue ?: 0.01
params.cutrun_qvalue = params.cutrun_qvalue ?: 0.05
params.min_cpg_depth = params.min_cpg_depth ?: 6
params.min_required_probe_fraction = params.min_required_probe_fraction ?: 0.80
params.sirt6_min_reciprocal_overlap = params.sirt6_min_reciprocal_overlap ?: 0.50

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
    output:
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.bam"), path("${meta.sample_id}.mapq${params.mapq}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.sample_id}.flagstat.txt"), emit: flagstat
    script:
    def readArgs = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    """
    set -euo pipefail
    bowtie2 -x ${params.bowtie2_index} ${readArgs} -p ${task.cpus} 2> ${meta.sample_id}.bowtie2.log | \
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
    output:
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.bam"), path("${meta.sample_id}.mapq${params.mapq}.bam.bai"), emit: bam
    script:
    def readArgs = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    """
    set -euo pipefail
    bowtie2 --very-sensitive --no-mixed --no-discordant --dovetail \
      -x ${params.bowtie2_index} ${readArgs} -p ${task.cpus} 2> ${meta.sample_id}.bowtie2.log | \
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

process ALIGN_WGBS {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/alignment/wgbs", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    // Paired-end methylation extraction requires mates to remain adjacent.
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.deduplicated.name.bam"), emit: bam
    tuple val(meta), path("${meta.sample_id}.wgbs_pair_filter.json"), emit: filter_stats
    script:
    def readArgs = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "${reads[0]}"
    def dedupMode = meta.paired ? '--paired' : '--single'
    def pairArg = meta.paired ? '--paired' : ''
    """
    set -euo pipefail
    mkdir -p bismark_${meta.sample_id}
    bismark --genome ${params.bismark_index} --bowtie2 --parallel 2 \
      --output_dir bismark_${meta.sample_id} ${readArgs}
    BAM=\$(find bismark_${meta.sample_id} -name '*_bismark_bt2*.bam' | head -n1)
    test -n "\${BAM}"
    deduplicate_bismark --bam ${dedupMode} --output_dir . \${BAM}
    DEDUP=\$(find . -maxdepth 1 -name '*.deduplicated.bam' | head -n1)
    test -n "\${DEDUP}"
    samtools sort -n -@ ${task.cpus} -o ${meta.sample_id}.deduplicated.name.pre-mapq.bam \${DEDUP}
    python3 ${projectDir}/bin/filter_paired_bam.py \
      --input ${meta.sample_id}.deduplicated.name.pre-mapq.bam \
      --output ${meta.sample_id}.mapq${params.mapq}.deduplicated.name.bam \
      --mapq ${params.mapq} ${pairArg} \
      --stats ${meta.sample_id}.wgbs_pair_filter.json
    samtools quickcheck -v ${meta.sample_id}.mapq${params.mapq}.deduplicated.name.bam
    rm -f ${meta.sample_id}.deduplicated.name.pre-mapq.bam
    """
    stub:
    """
    touch ${meta.sample_id}.mapq${params.mapq}.deduplicated.name.bam
    echo '{"mapq_threshold":${params.mapq},"paired":${meta.paired}}' > ${meta.sample_id}.wgbs_pair_filter.json
    """
}

process EXTRACT_METHYLATION {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/methylation", mode: 'copy'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("${meta.sample_id}.bismark.cov.gz"), emit: cov
    script:
    def pairedArg = meta.paired ? '--paired-end' : '--single-end'
    """
    set -euo pipefail
    bismark_methylation_extractor ${pairedArg} --comprehensive --gzip --bedGraph ${bam}
    COV=\$(find . -name '*.bismark.cov.gz' | head -n1)
    test -n "\${COV}"
    mv \${COV} ${meta.sample_id}.bismark.cov.gz
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
    echo 'sample_id,matched_probes,required_probes,fraction,min_cpg_depth' > pace_probe_qc.csv
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
    path 'manifest.json'
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
    echo '{}' > manifest.json
    """
}

workflow {
    samples = Channel.fromPath(params.samplesheet, checkIfExists: true).splitCsv(header: true).map { row ->
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

    all_reads = samples.mix(controls).mix(cutrun)
    trimmed = FASTP(all_reads).trimmed

    chip_and_input_trim = trimmed.filter { meta, reads -> meta.assay == 'CHIP' || meta.assay == 'INPUT' }
    wgbs_trim = trimmed.filter { meta, reads -> meta.assay == 'WGBS' }
    cutrun_trim = trimmed.filter { meta, reads -> meta.assay == 'CUTRUN' }

    chip_all_bam = ALIGN_CHIP(chip_and_input_trim).bam
    chip_bam = chip_all_bam.filter { meta, bam, bai -> meta.assay == 'CHIP' }
    input_bam = chip_all_bam.filter { meta, bam, bai -> meta.assay == 'INPUT' }
    control_bams = input_bam.map { meta, bam, bai -> bam }.collect()
    peaks = CALL_PEAKS(chip_bam, control_bams).peaks

    cutrun_all_bam = ALIGN_CUTRUN(cutrun_trim).bam
    cutrun_target_bam = cutrun_all_bam.filter { meta, bam, bai -> meta.target == 'SIRT6' }
    cutrun_igg_bam = cutrun_all_bam.filter { meta, bam, bai -> meta.target == 'IgG' }
    cutrun_controls = cutrun_igg_bam.map { meta, bam, bai -> bam }.collect()
    cutrun_peaks = CALL_CUTRUN_PEAKS(cutrun_target_bam, cutrun_controls).peaks
    sirt6_loci = BUILD_SIRT6_LOCI(cutrun_peaks.map { meta, peak -> peak }.collect()).loci

    wgbs_bam = ALIGN_WGBS(wgbs_trim).bam
    cov = EXTRACT_METHYLATION(wgbs_bam).cov
    pace_result = BUILD_PACE_MATRIX(cov.map { meta, path -> path }.collect())

    occupancy = CHIP_OCCUPANCY(
        peaks.map { meta, path -> path }.collect(),
        chip_bam.map { meta, bam, bai -> bam }.collect(),
        chip_bam.map { meta, bam, bai -> bai }.collect(),
        sirt6_loci
    ).occupancy

    corr = CORRELATE(pace_result.scores, occupancy).correlations
    MANIFEST(corr, pace_result.model_metadata)
}
