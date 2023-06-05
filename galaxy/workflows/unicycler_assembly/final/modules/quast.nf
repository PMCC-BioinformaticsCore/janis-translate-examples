nextflow.enable.dsl=2

process QUAST {
    debug true
    container "quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7"
    publishDir "${params.outdir}/quast"

    input:
    path assembly

    output:
    path "outputdir/report.html", emit: outReportHtml
    path "outputdir/report.pdf", emit: outReportPdf
    path "outputdir/report.tsv", emit: outReportTabular

    script:
    """
    quast \
    --ambiguity-score 0.99 \
    --ambiguity-usage "one" \
    --contig-thresholds "0,1000" \
    --extensive-mis-size 1000 \
    --gene-thresholds "0,300,1500,3000" \
    --k-mer-size 101 \
    --min-alignment 65 \
    --min-contig 500 \
    --min-identity 95.0 \
    --scaffold-gap-max-size 1000 \
    --threads 1 \
    --unaligned-part-size 500 \
    -o "outputdir" \
    ${assembly}
    """

}
