nextflow.enable.dsl=2

process RSEQC_GENE_BODY_COVERAGE {
    
    container "quay.io/biocontainers/rseqc:2.6.4--py27hf8a1672_2"
    publishDir "${params.outdir}/rseqc_gene_body_coverage"

    input:
    path batch_mode_input
    path refgene

    output:
    path "${batch_mode_input.simpleName}.geneBodyCoverage.curves.pdf", emit: outputcurvespdf
    path "${batch_mode_input.simpleName}.geneBodyCoverage.txt", emit: outputtxt

    script:
    """
    geneBody_coverage.py \
    -r ${refgene} \
    -i ${batch_mode_input} \
    -o ${batch_mode_input.simpleName} \
    """

}
