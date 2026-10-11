process SPECIFICITY_PLOTS {
    tag { specificity_csv.baseName }

    conda "${projectDir}/envs/specificity_plots.yml"

    publishDir "${params.outdir}/specificity_plots", mode: 'copy'

    input:
    path specificity_csv

    output:
    path 'specificity_groups_heatmap.svg', emit: heatmap

    script:
    """
    Rscript "${projectDir}/scripts/specificity_plots.R" \
        --input "${specificity_csv}" \
        --output specificity_groups_heatmap.svg
    """
}