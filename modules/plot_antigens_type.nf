
process PLOT_ANTIGENS_TYPE {

    tag "${specificity_csv.simpleName}"

    publishDir "${params.outdir}/antigens_type",
        mode: 'copy',
        overwrite: true

    input:
    path specificity_csv

    output:
    path "specificity_groups_by_type_heatmap.svg", emit: heatmap
    path "merged_summed_genes.csv", emit: merged

    script:
    """
    Rscript "${projectDir}/scripts/plot_antigens_type.R" \
        --input "${specificity_csv}" \
        --output specificity_groups_by_type_heatmap.svg \
        --merged_output merged_summed_genes.csv
    """
}
