nextflow.enable.dsl = 2

process FASTQC {
    
    container "quay.io/biocontainers/fastqc:0.11.8--2"
    publishDir "${params.outdir}/fastqc"

    input:
    path input_file
    path contaminants, stageAs: 'contaminants/*'
    path limits, stageAs: 'limits/*'
    val format_string
    val outdir

    output:
    path "${outdir}/*_fastqc.html", emit: out_html_file
    path "${outdir}/*_fastqc.zip", emit: out_text_file

    script:
    def contaminants = contaminants.simpleName != params.NULL_VALUE ? "--contaminants ${contaminants}" : ""
    def limits = limits.simpleName != params.NULL_VALUE ? "--limits ${limits}" : ""
    """
    mkdir ${outdir}
    fastqc \
    ${contaminants} \
    ${limits} \
    --outdir ${outdir} \
    -f ${format_string} \
    ${input_file} \
    """

}
