

ch_bam = Channel.fromPath( params.bam ).toList()

workflow {
    SAMTOOLS_FLAGSTAT(ch_bam)
}

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    publishDir "./outputs"

    input:
    path bam

    output:
    path "${bam[0]}.flagstat", emit: flagstats

    script:
    def bam = bam[0]
    """
    /usr/local/bin/samtools flagstat \
    ${bam} \
    > ${bam}.flagstat
    """

}
