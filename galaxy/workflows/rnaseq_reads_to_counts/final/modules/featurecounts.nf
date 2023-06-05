nextflow.enable.dsl=2

process FEATURECOUNTS {
    
    container "quay.io/biocontainers/janis-translate-featurecounts-2.0.1"
    publishDir "${params.outdir}/featurecounts"

    input:
    path alignment

    output:
    path "${alignment.simpleName}.txt", emit: output_short
    path "${alignment.simpleName}.txt.summary", emit: output_summary

    script:
    """
    featureCounts \
    -a /usr/local/annotation/mm10_RefSeq_exon.txt \
    -F "SAF" \
    -o ${alignment.simpleName}.txt \
    ${alignment} \
    """

}
