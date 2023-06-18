nextflow.enable.dsl = 2

process RSEQC_READ_DISTRIBUTION {
    
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_read_distribution"

    input:
    path input_bam
    path refgene

    output:
    path "${input_bam.simpleName}.read_distribution.txt", emit: outputTextfile

    script:
    """
    read_distribution.py \
    -i ${input_bam} \
    -r ${refgene} \
    > ${input_bam.simpleName}.read_distribution.txt \
    """

}
