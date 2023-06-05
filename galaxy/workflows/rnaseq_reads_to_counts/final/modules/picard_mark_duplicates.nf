nextflow.enable.dsl=2

process PICARD_MARK_DUPLICATES {
    
    container "quay.io/biocontainers/picard:2.18.2--py36_0"
    publishDir "${params.outdir}/picard_mark_duplicates"

    input:
    path input_file

    output:
    path "${input_file.simpleName}.markdup.bam", emit: outFile2
    path "${input_file.simpleName}.metrics.txt", emit: out_metrics_file

    script:
    """
    picard MarkDuplicates \
    INPUT=${input_file} \
    OUTPUT=${input_file.simpleName}.markdup.bam \
    METRICS_FILE=${input_file.simpleName}.metrics.txt \
    """

}
