nextflow.enable.dsl=2

process MARK_DUPLICATES_AND_SORT {
    
    container "mgibio/mark_duplicates-cwl:1.0.1"
    publishDir "${params.outdir}/mark_duplicates_and_sort"
    cpus "${params.mark_duplicates_and_sort.cpus}"
    memory "${params.mark_duplicates_and_sort.memory}"

    input:
    path script, stageAs: 'markduplicates_helper.sh'
    path bam
    val input_sort_order
    val output_name

    output:
    path "${bam.simpleName}.mark_dups_metrics.txt", emit: metrics_file
    tuple path(output_name), path("*.bai"), emit: sorted_bam

    script:
    def input_sort_order = input_sort_order != params.NULL_VALUE ? input_sort_order : "queryname"
    def output_name = output_name != params.NULL_VALUE ? output_name : "MarkedSorted.bam"
    """
    /bin/bash markduplicates_helper.sh \
    ${bam} \
    8 \
    ${output_name} \
    ${bam.simpleName}.mark_dups_metrics.txt \
    ${input_sort_order} \
    """

}
