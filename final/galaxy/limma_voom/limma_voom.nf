nextflow.enable.dsl = 2

workflow {
    LIMMA_VOOM(
        params.annotation_file,     // anno_geneanno
        params.contrast_file,       // cont_cinfo
        params.matrix_file,         // input_counts
        params.factor_file,         // input_fact_finfo
        params.limma_voom_script,   // limma_voom_script
        params.html_path,           // out_report1
        params.output_path,         // out_report_files_path
    )
}

process LIMMA_VOOM {
    
    container "quay.io/grace_hall1/limma-voom:3.50.1"
    publishDir "./outputs"

    input:
    path anno_geneanno
    path cont_cinfo
    path input_counts
    path input_fact_finfo
    path limma_voom_script
    val out_report1
    val out_report_files_path

    output:
    path "outReport.html", emit: outReport2

    script:
    """
    Rscript \
    ${limma_voom_script} \
    -C ${cont_cinfo} \
    -R ${out_report1} \
    -a ${anno_geneanno} \
    -f ${input_fact_finfo} \
    -m ${input_counts} \
    -G 10 \
    -P "i" \
    -c 1 \
    -d "BH" \
    -l 0 \
    -n "TMM" \
    -o ${out_report_files_path} \
    -p 0.05 \
    -s 0 \
    -t 3 \
    -z 0
    """

}
