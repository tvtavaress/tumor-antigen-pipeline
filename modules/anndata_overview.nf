process ANNDATA_OVERVIEW  {

    tag { h5ad_file.name }

    publishDir "${params.outdir}/anndata_overview", mode: 'copy'

    input:
    path h5ad_file

    output:
    path 'anndata_overview.h5ad', emit: h5ad
    path 'anndata_overview_report.html', emit: report

    script:
    """
    python "${projectDir}/scripts/anndata_overview.py" \
        --input "${h5ad_file}" \
        --output anndata_overview.h5ad \
        --report anndata_overview_report.html
    """
}