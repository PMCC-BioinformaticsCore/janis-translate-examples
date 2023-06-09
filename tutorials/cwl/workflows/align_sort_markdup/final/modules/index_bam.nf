nextflow.enable.dsl = 2

process INDEX_BAM {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    publishDir "${params.outdir}/index_bam"
    memory "${params.index_bam.memory}"

    input:
    path bam

    output:
    tuple path("${bam}"), path("${bam}.bai"), emit: indexed_bam

    script:
    """
    /usr/local/bin/samtools \
    index \
    ${bam} \
    ${bam}.bai \
     &&  \
    cp \
    ${bam}.bai \
    ${bam.simpleName}.bai
    """

}
