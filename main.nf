
nextflow.enable.dsl=2

include { ANNDATA_OVERVIEW } from './modules/anndata_overview'
include { SPECIFICITY } from './modules/specificity'
include { SPECIFICITY_PLOTS } from './modules/specificity_plots'


workflow {

    input_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )

    ANNDATA_OVERVIEW(input_ch)

    if (params.specificity == 'hpa') {

        SPECIFICITY(input_ch)

        SPECIFICITY_PLOTS(
            SPECIFICITY.out.summary
        )

    } else {
        error "Unsupported specificity method: ${params.specificity}"
    }
}
