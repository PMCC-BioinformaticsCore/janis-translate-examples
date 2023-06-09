nextflow.enable.dsl = 2

process NANOPLOT {
    debug true
    container "quay.io/biocontainers/nanoplot:1.28.2--py_0"
    publishDir "${params.outdir}/nanoplot"

    input:
    path fastq

    output:
    path "LogTransformed_HistogramReadlength.*", emit: outLogReadLength
    path "NanoStats.txt", emit: outNanostats
    path "HistogramReadlength.*", emit: outReadLength
    path "NanoPlot-report.html", emit: outputHtml

    script:
    """
    NanoPlot \
    --color "aliceblue" \
    --format "png" \
    --plots "kde" \
    --readtype "1D" \
    --threads 4 \
    -o "." \
    --fastq ${fastq}
    """

}
