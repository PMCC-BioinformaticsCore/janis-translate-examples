nextflow.enable.dsl=2

ch_bam = Channel.fromPath( params.bam )

workflow {
    SAMTOOLS_FLAGSTAT(
        ch_bam,             // input1 
        params.threads      // addthreads
    )
    SAMTOOLS_FLAGSTAT.out.output1.view()
}

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.13--h8c37831_0"

    input:
    path input1
    val addthreads

    output:
    path "output1.txt", emit: output1

    script:
    """
    samtools flagstat \
    --output-fmt "txt" \
    -@ ${addthreads} \
    ${input1} \
    > output1.txt
    """

}
