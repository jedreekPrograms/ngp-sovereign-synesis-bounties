nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: 'resources/samplesheet.csv'
params.controls_sheet = params.controls_sheet ?: 'resources/chip_inputs.csv'
params.outdir = params.outdir ?: 'results'
params.mapq = params.mapq ?: 30
params.peak_qvalue = params.peak_qvalue ?: 0.01
params.min_cpg_depth = params.min_cpg_depth ?: 6

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
    bowtie2 -x ${params.bowtie2_index} ${readArgs} -p ${task.cpus} 2> ${meta.sample_id}.bowtie2.log | \
      samtools view -b -q ${params.mapq} - | samtools sort -@ ${task.cpus} -o ${meta.sample_id}.mapq${params.mapq}.bam
    samtools index ${meta.sample_id}.mapq${params.mapq}.bam
    samtools flagstat ${meta.sample_id}.mapq${params.mapq}.bam > ${meta.sample_id}.flagstat.txt
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
    test -f ${control}
    macs3 callpeak -t ${bam} -c ${control} -f ${macsFormat} -g hs -n ${meta.sample_id} -q ${params.peak_qvalue} --keep-dup auto -B
    """
    stub:
    """
    printf 'chr1\t100\t200\t${meta.sample_id}\t100\t.\t1\t1\t1\t50\n' > ${meta.sample_id}_peaks.narrowPeak
    printf 'chr1\t100\t200\t1\n' > ${meta.sample_id}_treat_pileup.bdg
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
    stub:
    """
    touch ${meta.sample_id}.mapq${params.mapq}.deduplicated.bam
    touch ${meta.sample_id}.mapq${params.mapq}.deduplicated.bam.bai
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
    bismark_methylation_extractor ${pairedArg} --comprehensive --gzip --bedGraph ${bam}
    COV=\$(find . -name '*.bismark.cov.gz' | head -n1)
    test -n "\${COV}"
    mv \${COV} ${meta.sample_id}.bismark.cov.gz
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
    python3 ${projectDir}/bin/build_manifest.py \
      --correlations ${correlations} \
      --model-metadata ${model_metadata} \
      --mapq ${params.mapq} \
      --peak-fdr ${params.peak_qvalue} \
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

    all_reads = samples.mix(controls)
    trimmed = FASTP(all_reads).trimmed
    chip_and_input_trim = trimmed.filter { meta, reads -> meta.assay == 'CHIP' || meta.assay == 'INPUT' }
    wgbs_trim = trimmed.filter { meta, reads -> meta.assay == 'WGBS' }

    chip_all_bam = ALIGN_CHIP(chip_and_input_trim).bam
    chip_bam = chip_all_bam.filter { meta, bam, bai -> meta.assay == 'CHIP' }
    input_bam = chip_all_bam.filter { meta, bam, bai -> meta.assay == 'INPUT' }
    control_bams = input_bam.map { meta, bam, bai -> bam }.collect()

    peaks = CALL_PEAKS(chip_bam, control_bams).peaks
    wgbs_bam = ALIGN_WGBS(wgbs_trim).bam
    cov = EXTRACT_METHYLATION(wgbs_bam).cov
    pace_result = BUILD_PACE_MATRIX(cov.map { meta, p -> p }.collect())
    occupancy = CHIP_OCCUPANCY(
        peaks.map { meta, p -> p }.collect(),
        chip_bam.map { meta, b, bai -> b }.collect(),
        chip_bam.map { meta, b, bai -> bai }.collect()
    ).occupancy
    corr = CORRELATE(pace_result.scores, occupancy).correlations
    MANIFEST(corr, pace_result.model_metadata)
}
