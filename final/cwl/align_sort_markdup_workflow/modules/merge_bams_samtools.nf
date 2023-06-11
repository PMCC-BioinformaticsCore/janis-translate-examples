nextflow.enable.dsl = 2

process MERGE_BAMS_SAMTOOLS {
    
    container "quay.io/biocontainers/samtools:0.1.19--hfb9b9cc_8"
    publishDir "${params.outdir}/merge_bams_samtools"
    cpus "${params.merge_bams_samtools.cpus}"
    memory "${params.merge_bams_samtools.memory}"

    input:
    path bams
    val name

    output:
    path "${name}.merged.bam", emit: merged_bam

    script:
    def bams_joined = bams.join(' ')
    """
    /usr/local/bin/samtools merge \
    -@ 4 \
    ${name}.merged.bam \
    ${bams_joined}
    """

}
