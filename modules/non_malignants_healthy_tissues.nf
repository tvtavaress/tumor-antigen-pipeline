process NON_MALIGNANTS_HEALTHY_TISSUES {

    tag 'HPA_association'
    
    publishDir "${params.outdir}/non_malignants_healthy_tissues", mode: 'copy'

    conda "${projectDir}/envs/analysis.yml"

    input:
    path cta_file

    output:
    path 'association.tsv', emit: association

    script:
    """
    python "${projectDir}/scripts/non_malignants_healthy_tissues.py" \
        --cta "${cta_file}" \
        --output association.tsv
    """
}