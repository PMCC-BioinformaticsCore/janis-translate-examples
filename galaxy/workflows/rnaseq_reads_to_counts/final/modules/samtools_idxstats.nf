nextflow.enable.dsl=2

process SAMTOOLS_IDXSTATS {
    
    container "quay.io/biocontainers/samtools:1.9--h10a08f8_12"
    publishDir "${params.outdir}/samtools_idxstats"

    input:
    path input_bam
    val addthreads

    output:
    path "${input_bam.simpleName}.idxstats.tabular", emit: outputTabular

    script:
    """
    samtools idxstats \
    -@ ${addthreads} \
    ${input_bam} \
    > ${input_bam.simpleName}.idxstats.tabular \
    """

}
