nextflow.enable.dsl=2

process COLLECTION_COLUMN_JOIN {
    
    container "quay.io/biocontainers/coreutils:8.25--1"
    publishDir "${params.outdir}/collection_column_join"

    input:
    path input_tabular
    path collection_column_join_script

    output:
    path "unknown_collection_pattern", emit: out_tabular_output

    script:
    def input_tabular_joined = input_tabular.join(' ')
    """
    sh \
    ${input_tabular_joined} \
    ${collection_column_join_script} \
    """

}
