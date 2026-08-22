nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: 'resources/samplesheet.csv'
params.outdir = params.outdir ?: 'results'
params.mapq = params.mapq ?: 30
params.peak_qvalue = params.peak_qvalue ?: 0.01

process FASTP {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/qc/fastp", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.sample_id}.trimmed*.fastq.gz"), emit: trimmed
    path "${meta.sample_id}.fastp.json", emit: json
    path "${meta.sample_id}.fastp.html", emit: html
    script:
    def ioArgs = reads.size() == 2 ? "-i ${reads[0]} -I ${reads[1]} -o ${meta.sample_id}.trimmed.R1.fastq.gz -O ${meta.sample_id}.trimmed.R2.fastq.gz" : "-i ${reads[0]} -o ${meta.sample_id}.trimmed.R1.fastq.gz"
    """
    fastp ${ioArgs} --json ${meta.sample_id}.fastp.json --html ${meta.sample_id}.fastp.html --thread ${task.cpus}
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
    bowtie2 -x ${params.bowtie2_index} ${readArgs} -p ${task.cpus} 2> ${meta.sample_id}.bowtie2.log | \
      samtools view -b -q ${params.mapq} - | samtools sort -@ ${task.cpus} -o ${meta.sample_id}.mapq${params.mapq}.bam
    samtools index ${meta.sample_id}.mapq${params.mapq}.bam
    samtools flagstat ${meta.sample_id}.mapq${params.mapq}.bam > ${meta.sample_id}.flagstat.txt
    """
}

process CALL_PEAKS {
    tag "${meta.sample_id}:${meta.mark}"
    publishDir "${params.outdir}/peaks", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.sample_id}_peaks.narrowPeak"), emit: peaks
    tuple val(meta), path("${meta.sample_id}_treat_pileup.bdg"), emit: pileup
    script:
    def macsFormat = meta.paired ? 'BAMPE' : 'BAM'
    """
    macs3 callpeak -t ${bam} -f ${macsFormat} -g hs -n ${meta.sample_id} -q ${params.peak_qvalue} --keep-dup auto -B
    """
}

process ALIGN_WGBS {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/alignment/wgbs", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.sample_id}.mapq${params.mapq}.deduplicated.bam"), path("${meta.sample_id}.mapq${params.mapq}.deduplicated.bam.bai"), emit: bam
    script:
    def readArgs = reads.size() == 2 ? "-1 ${reads[0]} -2 ${reads[1]}" : "${reads[0]}"
    """
    mkdir -p bismark_${meta.sample_id}
    bismark --genome ${params.bismark_index} --bowtie2 --parallel 2 --output_dir bismark_${meta.sample_id} ${readArgs}
    BAM=\$(find bismark_${meta.sample_id} -name '*_bismark_bt2*.bam' | head -n1)
    test -n "\${BAM}"
    deduplicate_bismark --bam ${meta.paired ? '--paired' : ''} --output_dir . \${BAM}
    DEDUP=\$(find . -maxdepth 1 -name '*.deduplicated.bam' | head -n1)
    test -n "\${DEDUP}"
    samtools view -b -q ${params.mapq} \${DEDUP} | samtools sort -@ ${task.cpus} -o ${meta.sample_id}.mapq${params.mapq}.deduplicated.bam
    samtools index ${meta.sample_id}.mapq${params.mapq}.deduplicated.bam
    """
}

process EXTRACT_METHYLATION {
    tag "${meta.sample_id}"
    label 'big_mem'
    publishDir "${params.outdir}/methylation", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.sample_id}.bismark.cov.gz"), emit: cov
    script:
    def pairedArg = meta.paired ? '--paired-end' : ''
    """
    bismark_methylation_extractor ${pairedArg} --gzip --bedGraph --cytosine_report --genome_folder ${params.bismark_index} ${bam}
    COV=\$(find . -name '*.bismark.cov.gz' | head -n1)
    test -n "\${COV}"
    mv \${COV} ${meta.sample_id}.bismark.cov.gz
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
    script:
    """
    Rscript ${projectDir}/bin/compute_dunedinpace.R --cov-dir . --scores dunedinpace_scores.csv --qc pace_probe_qc.csv --min-probe-fraction ${params.min_required_probe_fraction}
    """
}

process CHIP_OCCUPANCY {
    tag 'consensus-occupancy'
    publishDir "${params.outdir}/chip", mode: 'copy'
    input:
    path peak_files
    path bam_files
    path bai_files
    output:
    path 'histone_occupancy.csv', emit: occupancy
    script:
    """
    python3 ${projectDir}/bin/compute_chip_occupancy.py --peaks-glob '*_peaks.narrowPeak' --bam-glob '*.bam' --output histone_occupancy.csv
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
    python3 ${projectDir}/bin/correlate.py --pace ${scores} --occupancy ${occupancy} --output correlations.json
    """
}

process MANIFEST {
    tag 'manifest'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path correlations
    output:
    path 'manifest.json'
    script:
    """
    python3 ${projectDir}/bin/build_manifest.py --correlations ${correlations} --mapq ${params.mapq} --peak-fdr ${params.peak_qvalue} --output manifest.json
    """
}

workflow {
    samples = Channel.fromPath(params.samplesheet, checkIfExists: true).splitCsv(header: true).map { row ->
        def meta = [sample_id: row.sample_id, condition: row.condition, replicate: row.replicate as Integer, assay: row.assay, mark: row.mark ?: 'NA', paired: row.fastq_2 ? true : false]
        def reads = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
        tuple(meta, reads)
    }
    chip_reads = samples.filter { meta, reads -> meta.assay == 'CHIP' }
    wgbs_reads = samples.filter { meta, reads -> meta.assay == 'WGBS' }
    chip_trim = FASTP(chip_reads).trimmed
    wgbs_trim = FASTP(wgbs_reads).trimmed
    chip_bam = ALIGN_CHIP(chip_trim).bam
    peaks = CALL_PEAKS(chip_bam).peaks
    wgbs_bam = ALIGN_WGBS(wgbs_trim).bam
    cov = EXTRACT_METHYLATION(wgbs_bam).cov
    pace = BUILD_PACE_MATRIX(cov.map { meta, p -> p }.collect()).scores
    occupancy = CHIP_OCCUPANCY(peaks.map { meta, p -> p }.collect(), chip_bam.map { meta, b, bai -> b }.collect(), chip_bam.map { meta, b, bai -> bai }.collect()).occupancy
    corr = CORRELATE(pace, occupancy).correlations
    MANIFEST(corr)
}
