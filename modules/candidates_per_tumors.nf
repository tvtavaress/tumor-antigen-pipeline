process CANDIDATES_PER_TUMORS {

    tag "candidates_per_tumors"

    publishDir "${params.outdir}/candidates_per_tumors", mode: 'copy'

    input:
    path h5ad
    val gene_name_mapping
    path candidate_genes

    output:
    path "candidate_gene_pct_expression_by_tumor.csv", emit: percentages
    path "candidate_gene_binary_matrix.csv", emit: binary_matrix
    path "histogram_all_genes.png", emit: histogram_all_genes
    path "histogram_candidate_genes.png", emit: histogram_candidates

    script:
    def mapping_arg = gene_name_mapping
        ? "--gene-name-mapping '${gene_name_mapping}'"
        : ""

    """
    python "${projectDir}/scripts/candidates_per_tumors.py" \
        --input "${h5ad}" \
        --candidate-genes "${candidate_genes}" \
        --output-dir "." \
        --celltype-column "cell_type" \
        --malignant-label "malignant cell" \
        --tumor-column "disease" \
        --candidate-gene-column "gene_name" \
        --expression-threshold 0.1 \
        --binary-threshold 1 \
        ${mapping_arg}
    """
}