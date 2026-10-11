process SET_ANTIGENS_TYPE {

    publishDir "${params.outdir}/antigens_type", mode: 'copy'

    input:
    path summary_csv

    output:
    path "classified_specificity_labeled.csv", emit: classified_labeled

    script:
    """
    python3 ${projectDir}/scripts/set_antigens_type.py \
        --input ${summary_csv}
    """
}