process SPECIFICITY {

    tag "${input_file.simpleName}"

    publishDir "${params.outdir}/specificity", mode: 'copy'

    input:
    path input_file

    output:
    path 'classified_specificity_*.csv', emit: classified
    path 'specificity_summary_*.csv', emit: summary
    path 'enriched_malignant_cell_*.csv', emit: enriched

    script:
    """
    python "${projectDir}/scripts/specificity.py" \
        --input "${input_file}" \
        --output-dir .
    """
}