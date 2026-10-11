process PLOT_SPECIFICITY {
    tag { specificity_csv.baseName }

    conda "${projectDir}/envs/specificity_plots.yml"

    publishDir "${params.outdir}/specificity", mode: 'copy'

    input:
    path specificity_csv

    output:
    path 'specificity_groups_heatmap.svg', emit: heatmap

    script:
    """
    Rscript "${projectDir}/scripts/plot_specificity.R" \
        --input "${specificity_csv}" \
        --output specificity_groups_heatmap.svg
    """
}