nextflow.enable.dsl=2

process MULTIQC {
    
    container "quay.io/biocontainers/multiqc:1.6--py35h24bf2e0_0"
    publishDir "${params.outdir}/multiqc"

    input:
    path multiqc_config
    path unknown1
    path unknown2
    path unknown3
    path unknown4
    path unknown5
    path unknown6
    path unknown7
    path unknown8
    path unknown9

    output:
    path "multiqc_report.html", emit: out_html_report
    path "multiqc_data/multiqc_*.txt", emit: out_stats

    script:
    """
    multiqc . \
    """

}
