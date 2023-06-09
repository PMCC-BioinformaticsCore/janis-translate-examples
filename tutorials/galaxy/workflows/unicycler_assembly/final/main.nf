nextflow.enable.dsl = 2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { UNICYCLER } from './modules/unicycler'
include { NANOPLOT } from './modules/nanoplot'
include { QUAST } from './modules/quast'
include { BUSCO } from './modules/busco'


// data which will be passed as channels
ch_in_long      = Channel.fromPath( params.in_long )
ch_in_short_r1  = Channel.fromPath( params.in_short_r1 )
ch_in_short_r2  = Channel.fromPath( params.in_short_r2 )


// data which will be passed as optional files
fastqc1_adapters         = file( params.fastqc1_adapters )
fastqc1_contaminants     = file( params.fastqc1_contaminants )
fastqc1_limits           = file( params.fastqc1_limits )
fastqc2_adapters         = file( params.fastqc2_adapters )
fastqc2_contaminants     = file( params.fastqc2_contaminants )
fastqc2_limits           = file( params.fastqc2_limits )
unicycler_contamination  = file( params.unicycler_contamination )
unicycler_start_genes    = file( params.unicycler_start_genes )


workflow {

    FASTQC1(
        ch_in_short_r1,
        fastqc1_adapters,
        fastqc1_contaminants,
        fastqc1_limits
    )

    FASTQC2(
        ch_in_short_r2,
        fastqc2_adapters,
        fastqc2_contaminants,
        fastqc2_limits
    )

    UNICYCLER(
        ch_in_short_r1,
        ch_in_short_r2,
        unicycler_contamination,
        ch_in_long,
        unicycler_start_genes
    )

    NANOPLOT(
        ch_in_long
    )

    QUAST(
        UNICYCLER.out.outAssembly
    )

    BUSCO(
        UNICYCLER.out.outAssembly
    )


}
