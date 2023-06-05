nextflow.enable.dsl=2

process RSEQC_INFER_EXPERIMENT {
    
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_infer_experiment"

    input:
    path input_bam
    path refgene

    output:
    path "${input_bam.simpleName}.infer_experiment.txt", emit: outputTextfile

    script:
    """
    infer_experiment.py \
    -i ${input_bam} \
    -r ${refgene} \
    > ${input_bam.simpleName}.infer_experiment.txt \
    """

}
