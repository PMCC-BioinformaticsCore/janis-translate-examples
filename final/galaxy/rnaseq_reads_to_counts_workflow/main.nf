nextflow.enable.dsl = 2

include { FASTQC as FASTQC1 } from './modules/fastqc'
include { CUTADAPT } from './modules/cutadapt'
include { FASTQC as FASTQC2 } from './modules/fastqc'
include { HISAT2 } from './modules/hisat2'
include { FEATURECOUNTS } from './modules/featurecounts'
include { PICARD_MARK_DUPLICATES } from './modules/picard_mark_duplicates'
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats'
include { RSEQC_GENE_BODY_COVERAGE } from './modules/rseqc_gene_body_coverage'
include { RSEQC_INFER_EXPERIMENT } from './modules/rseqc_infer_experiment'
include { RSEQC_READ_DISTRIBUTION } from './modules/rseqc_read_distribution'
include { COLLECTION_COLUMN_JOIN } from './modules/collection_column_join'
include { MULTIQC } from './modules/multiqc'


// data which will be passed as channels
ch_in_input_fastqs_collection  = Channel.fromPath( params.in_input_fastqs_collection ).toList()
ch_hisat2_index                = Channel.fromPath( params.hisat2_index_path ).toList()

// data which will be passed as variables
collection_column_join_script  = file( params.collection_column_join_script )
fastqc1_contaminants           = file( params.fastqc1_contaminants )
fastqc1_limits                 = file( params.fastqc1_limits )
fastqc2_contaminants           = file( params.fastqc2_contaminants )
fastqc2_limits                 = file( params.fastqc2_limits )
in_input_reference_gene_bed    = file( params.in_input_reference_gene_bed )
multiqc_config                 = file( params.multiqc_config )


workflow {

    FASTQC1(
        ch_in_input_fastqs_collection.flatten(),  // input_file
        fastqc1_contaminants,                     // contaminants
        fastqc1_limits,                           // limits
        params.fastqc1_format_string,             // format_string
        params.fastqc1_outdir                     // outdir
    )

    CUTADAPT(
        ch_in_input_fastqs_collection.flatten(),  // library_input_1
        params.cutadapt_adapter
    )

    FASTQC2(
        CUTADAPT.out.out12,            // input_file
        fastqc2_contaminants,          // contaminants
        fastqc2_limits,                // limits
        params.fastqc2_format_string,  // format_string
        params.fastqc2_outdir          // outdir
    )

    HISAT2(
        CUTADAPT.out.out12,        // library_input_1
        ch_hisat2_index,           // index_path
    )

    FEATURECOUNTS(
        HISAT2.out.output_alignments  // alignment
    )

    PICARD_MARK_DUPLICATES(
        HISAT2.out.output_alignments  // input_file
    )

    SAMTOOLS_IDXSTATS(
        HISAT2.out.output_alignments,        // input_bam
        params.samtools_idxstats_addthreads  // addthreads
    )

    RSEQC_GENE_BODY_COVERAGE(
        HISAT2.out.output_alignments,             // batch_mode_input
        in_input_reference_gene_bed,              // refgene
    )

    RSEQC_INFER_EXPERIMENT(
        HISAT2.out.output_alignments,  // input_bam
        in_input_reference_gene_bed    // refgene
    )

    RSEQC_READ_DISTRIBUTION(
        HISAT2.out.output_alignments,  // input_bam
        in_input_reference_gene_bed    // refgene
    )

    // COLLECTION_COLUMN_JOIN(
    //     FEATURECOUNTS.out.output_short.toList(),  // input_tabular
    //     collection_column_join_script             // collection_column_join_script
    // )

    MULTIQC(
        multiqc_config,                               // multiqc_config
        FASTQC2.out.out_text_file.toList(),                    // unknown1
        CUTADAPT.out.out_report.toList(),                      // unknown2
        RSEQC_INFER_EXPERIMENT.out.outputTextfile.toList(),    // unknown3
        PICARD_MARK_DUPLICATES.out.out_metrics_file.toList(),  // unknown4
        SAMTOOLS_IDXSTATS.out.outputTabular.toList(),          // unknown5
        RSEQC_GENE_BODY_COVERAGE.out.outputtxt.toList(),       // unknown6
        RSEQC_READ_DISTRIBUTION.out.outputTextfile.toList(),   // unknown7
        FEATURECOUNTS.out.output_summary.toList(),             // unknown8
        HISAT2.out.out_summary_file.toList()                   // unknown9
    )


}
