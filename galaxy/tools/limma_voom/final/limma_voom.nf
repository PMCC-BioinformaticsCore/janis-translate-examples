nextflow.enable.dsl=2

workflow {
    LIMMA_VOOM(
        params.limma_voom_script,   // limma_voom_script
        params.annotation_file,     // option_a
        params.contrast_file,       // option_c
        params.factor_file,         // option_f
        params.matrix_file,         // option_m
        params.html_path,           // option_r
        params.output_path,         // option_o
    )
}

nextflow.enable.dsl=2

process LIMMA_VOOM {
    
    container "quay.io/biocontainers/janis-translate-limma-voom-3.34.9.9"
    publishDir "./outputs"

    input:
    path limma_voom_script
    path option_a
    path option_c
    path option_f
    path option_m
    val option_r
    val option_o

    output:
    path "outReport.html", emit: outReport

    script:
    """
    Rscript \
    ${limma_voom_script} \
    -C ${option_c} \
    -R ${option_r} \
    -a ${option_a} \
    -f ${option_f} \
    -m ${option_m} \
    -G 10 \
    -P "i" \
    -c 1 \
    -d "BH" \
    -l 0 \
    -n "TMM" \
    -o ${option_o} \
    -p 0.05 \
    -s 0 \
    -t 3 \
    -z 0 \
    """

}

// process LIMMA_VOOM {
    
//     container "quay.io/biocontainers/janis-translate-limma-voom-3.34.9.9"
//     publishDir "./outputs"

//     input:
//     path limma_voom_script
//     path option_a
//     path option_c
//     path option_f
//     path option_m
//     val option_r
//     val option_o

//     output:
//     path "outReport.html", emit: outReport

//     script:
//     """
//     Rscript \
//     ${limma_voom_script} \
//     -C ${option_c} \
//     -R ${option_r} \
//     -a ${option_a} \
//     -f ${option_f} \
//     -m ${option_m} \
//     -G 10 \
//     -c 1 \
//     -d "BH" \
//     -l 0 \
//     -n "TMM" \
//     -o ${option_o} \
//     -p 0.05 \
//     -s 0 \
//     -t 3 \
//     -z 0 \
//     """

// }
