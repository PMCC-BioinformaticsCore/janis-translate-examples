nextflow.enable.dsl=2

process UNICYCLER {
    debug true
    container "quay.io/biocontainers/unicycler:0.4.8--py39h98c8e45_5"
    publishDir "${params.outdir}/unicycler"

    input:
    path option11
    path option12
    path contamination, stageAs: 'contamination/*'
    path option_l, stageAs: 'option_l/*'
    path start_genes, stageAs: 'start_genes/*'

    output:
    path "assembly.fasta", emit: outAssembly
    path "assembly.gfa", emit: outAssemblyGraph

    script:
    def contamination = contamination.simpleName != params.NULL ? "--contamination ${contamination}" : ""
    def option_l = option_l.simpleName != params.NULL ? "-l ${option_l}" : ""
    def start_genes = start_genes.simpleName != params.NULL ? "--start_genes ${start_genes}" : ""
    """
    unicycler \
    ${contamination} \
    ${start_genes} \
    -1 ${option11} \
    -2 ${option12} \
    ${option_l} \
    --depth_filter 0.25 \
    --kmer_count 10 \
    --linear_seqs 0 \
    --max_kmer_frac 0.95 \
    --min_component_size 1000 \
    --min_dead_end_size 1000 \
    --min_fasta_length 500 \
    --min_kmer_frac 0.2 \
    --min_polish_size 1000 \
    --mode "normal" \
    --scores "3,-6,-5,-2" \
    --start_gene_cov 95.0 \
    --start_gene_id 90.0 \
    --verbosity 3 \
    -o "./" \
    -t 4 \
    """

}
