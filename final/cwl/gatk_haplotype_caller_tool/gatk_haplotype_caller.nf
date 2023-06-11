
nextflow.enable.dsl = 2

ch_bam = Channel.fromPath( params.bam ).toList()
ch_dbsnp = Channel.fromPath( params.dbsnp ).toList()
ch_reference = Channel.fromPath( params.reference ).toList()

workflow {

    GATK_HAPLOTYPE_CALLER(
        ch_bam,
        ch_reference,
        ch_dbsnp,
        params.intervals,
        params.gvcf_gq_bands,
        params.emit_reference_confidence,
        params.contamination_fraction,
        params.max_alternate_alleles,
        params.ploidy,
        params.read_filter,
    )

}

process GATK_HAPLOTYPE_CALLER {
    
    container "broadinstitute/gatk:4.1.8.1"

    input:
    path bam
    path reference
    path dbsnp_vcf, stageAs: 'dbsnp_vcf/*'
    val intervals
    val gvcf_gq_bands
    val emit_reference_confidence
    val contamination_fraction
    val max_alternate_alleles
    val ploidy
    val read_filter

    output:
    tuple path("output.g.vcf.gz"), path("*.tbi"), emit: gvcf

    script:
    def bam = bam[0]
    def dbsnp_vcf = dbsnp_vcf[0] != null ? "--dbsnp ${dbsnp_vcf[0]}" : ""
    def reference = reference[0]
    def gvcf_gq_bands_joined = gvcf_gq_bands != params.NULL_VALUE ? "-GQB " + gvcf_gq_bands.join(' ') : ""
    def intervals_joined = intervals.join(' ')
    def contamination_fraction = contamination_fraction != params.NULL_VALUE ? "-contamination ${contamination_fraction}" : ""
    def max_alternate_alleles = max_alternate_alleles != params.NULL_VALUE ? "--max_alternate_alleles ${max_alternate_alleles}" : ""
    def ploidy = ploidy != params.NULL_VALUE ? "-ploidy ${ploidy}" : ""
    def read_filter = read_filter != params.NULL_VALUE ? "--read_filter ${read_filter}" : ""
    """
    /gatk/gatk --java-options -Xmx16g HaplotypeCaller \
    -R ${reference} \
    -I ${bam} \
    -ERC ${emit_reference_confidence} \
    ${gvcf_gq_bands_joined} \
    -L ${intervals_joined} \
    ${dbsnp_vcf} \
    ${contamination_fraction} \
    ${max_alternate_alleles} \
    ${ploidy} \
    ${read_filter} \
    -O "output.g.vcf.gz" \
    """

}
