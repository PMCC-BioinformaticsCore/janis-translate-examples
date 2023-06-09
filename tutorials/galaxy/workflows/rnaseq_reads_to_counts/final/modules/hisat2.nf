nextflow.enable.dsl = 2

process HISAT2 {
    
    container "quay.io/biocontainers/janis-translate-hisat2-2.2.1"
    publishDir "${params.outdir}/hisat2"

    input:
    path library_input_1
    path index_path

    output:
    path "${library_input_1.simpleName}.alignment_summary.txt", emit: out_summary_file
    path "${library_input_1.simpleName}.bam", emit: output_alignments

    script:
    """
    hisat2 \
    -U ${library_input_1} \
    -x ${index_path[0].simpleName} \
    --summary-file ${library_input_1.simpleName}.alignment_summary.txt \
    -S out.sam

    samtools view out.sam -o out.bam
    samtools sort out.bam -o sorted.bam
    samtools index sorted.bam -o sorted.bam.bai
    mv sorted.bam ${library_input_1.simpleName}.bam 
    mv sorted.bam.bai ${library_input_1.simpleName}.bam.bai
    """

}
