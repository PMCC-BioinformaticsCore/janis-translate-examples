nextflow.enable.dsl=2

process FASTQC {
    debug true
    container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
    publishDir "${params.outdir}/fastqc"

    input:
    path input_file
    path adapters, stageAs: 'adapters/*'
    path contaminants, stageAs: 'contaminants/*'
    path limits, stageAs: 'limits/*'

    output:
    path "${input_file.simpleName}_fastqc.html", emit: outHtmlFile

    script:
    def adapters = adapters.simpleName != params.NULL ? "--adapters ${adapters}" : ""
    def contaminants = contaminants.simpleName != params.NULL ? "--contaminants ${contaminants}" : ""
    def limits = limits.simpleName != params.NULL ? "--limits ${limits}" : ""
    """
    fastqc \
    ${adapters} \
    ${contaminants} \
    ${limits} \
    --kmers 7 \
    --threads 2 \
    ${input_file} \
    """

}
