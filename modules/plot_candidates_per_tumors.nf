process PLOT_CANDIDATES_PER_TUMORS {

    tag "candidate_gene_binary_matrix"

    publishDir "${params.outdir}/candidates_per_tumors", mode: 'copy'

    input:
    path binary_matrix

    output:
    path "upset_candidate_genes.svg", emit: plot
    path "genes_by_tumor_combination_summary.csv", emit: summary

    script:
    """
    Rscript "${projectDir}/scripts/plot_candidates_per_tumors.R" \
        --input "${binary_matrix}" \
        --output_dir "."
    """
}