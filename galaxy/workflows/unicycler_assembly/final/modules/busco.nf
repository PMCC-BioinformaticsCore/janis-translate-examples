nextflow.enable.dsl=2

process BUSCO {
    debug true
    container "quay.io/biocontainers/busco:5.3.2--pyhdfd78af_0"
    publishDir "${params.outdir}/busco"

    input:
    path in_file

    output:
    path "busco_galaxy/run_*/missing_busco_list.tsv", emit: outBuscoMissing
    path "busco_galaxy/run_*/short_summary.txt", emit: outBuscoSum
    path "busco_galaxy/run_*/full_table.tsv", emit: outBuscoTable

    script:
    """
    busco \
    --in ${in_file} \
    --augustus_species "local" \
    --cpu 4 \
    --evalue 0.001 \
    --limit 3 \
    --lineage_dataset "acidobacteria_odb10" \
    --mode "geno" \
    --out "busco_galaxy" \
    """

}
