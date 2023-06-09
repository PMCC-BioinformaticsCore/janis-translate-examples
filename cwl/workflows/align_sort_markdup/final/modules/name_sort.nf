nextflow.enable.dsl = 2

process NAME_SORT {
    
    container "mgibio/sambamba-cwl:0.6.4"
    publishDir "${params.outdir}/name_sort"
    cpus "${params.name_sort.cpus}"
    memory "${params.name_sort.memory}"

    input:
    path bam

    output:
    path "${bam.simpleName}.NameSorted.bam", emit: name_sorted_bam

    script:
    """
    /usr/bin/sambamba sort \
    -t \
    8 \
    -m \
    22G \
    -n \
    -o \
    ${bam.simpleName}.NameSorted.bam \
    ${bam}
    """

}
