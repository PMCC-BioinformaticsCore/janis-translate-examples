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

    output:
    path "refAlign.bam", emit: aligned_bam

    script:
    def reference = reference[0]
    """
    /bin/bash /usr/bin/alignment_helper.sh \
    ${bam} \
    "${readgroup}" \
    ${reference} \
    8 \
    > refAlign.bam \
    """

}
