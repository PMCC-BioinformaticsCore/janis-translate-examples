nextflow.enable.dsl=2

process ALIGN_AND_TAG {
    
    container "mgibio/alignment_helper-cwl:1.0.0"
    publishDir "${params.outdir}/align_and_tag"
    cpus "${params.align_and_tag.cpus}"
    memory "${params.align_and_tag.memory}"

    input:
    path reference
    path bam
    val readgroup
    val dummy

    output:
    path "${bam.simpleName}_refAlign.bam", emit: aligned_bam

    script:
    def reference = reference[0]
    def dummy = dummy != params.NULL_VALUE ? dummy : ""
    """
    /bin/bash /usr/bin/alignment_helper.sh \
    ${bam} \
    "${readgroup}" \
    ${reference} \
    ${dummy} \
    8 \
    > ${bam.simpleName}_refAlign.bam \
    """

}
