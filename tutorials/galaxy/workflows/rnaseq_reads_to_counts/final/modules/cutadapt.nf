nextflow.enable.dsl = 2

process CUTADAPT {
    
    container "quay.io/biocontainers/cutadapt:1.16--py35_2"
    publishDir "${params.outdir}/cutadapt"

    input:
    path library_input_1
    val adapter

    output:
    path "${library_input_1.simpleName}.cutadapt.fastq.gz", emit: out12
    path "${library_input_1.simpleName}.cutadapt.txt", emit: out_report

    script:
    """
    cutadapt \
    -a ${adapter} \
    -o ${library_input_1.simpleName}.cutadapt.fastq.gz \
    ${library_input_1} \
    > ${library_input_1.simpleName}.cutadapt.txt
    """

}

