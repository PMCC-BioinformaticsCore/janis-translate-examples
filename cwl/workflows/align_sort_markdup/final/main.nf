nextflow.enable.dsl=2

include { ALIGN } from './subworkflows/align'
include { INDEX_BAM } from './modules/index_bam'
include { MARK_DUPLICATES_AND_SORT } from './modules/mark_duplicates_and_sort'
include { MERGE_BAMS_SAMTOOLS as MERGE } from './modules/merge_bams_samtools'
include { NAME_SORT } from './modules/name_sort'


// data which will be passed as channels
ch_bams        = Channel.fromPath( params.bams ).toList()
ch_reference   = Channel.fromPath( params.reference ).toList()
ch_readgroups  = Channel.of( params.readgroups ).toList()


workflow {

    ALIGN(
        ch_bams.flatten(),       // bam
        ch_reference,            // reference
        ch_readgroups.flatten()  // readgroup
    )

    MERGE(
        ALIGN.out.tagged_bam.collect(),  // bams
        params.final_name               // name
    )

    NAME_SORT(
        MERGE.out.merged_bam  // bam
    )
    
    MARK_DUPLICATES_AND_SORT(
        params.mark_duplicates_and_sort.script,  // script
        NAME_SORT.out.name_sorted_bam,           // bam
        params.NULL_VALUE,                       // input_sort_order
        params.final_name                        // output_name
    )
    
    INDEX_BAM(
        MARK_DUPLICATES_AND_SORT.out.sorted_bam.map{ tuple -> tuple[0] }  // bam
    )

}
