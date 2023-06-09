nextflow.enable.dsl = 2

include { ALIGN_AND_TAG } from '../modules/align_and_tag'


workflow ALIGN {

    take:
    ch_bam
    ch_reference
    ch_readgroup

    main:
    ALIGN_AND_TAG(
        ch_reference[0, 1, 2, 3, 6, 7],  // reference
        ch_bam,                          // bam
        ch_readgroup                     // readgroup
    )

    emit:
    tagged_bam = ALIGN_AND_TAG.out.aligned_bam

}
