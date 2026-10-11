
nextflow.enable.dsl=2

include { ANNDATA_OVERVIEW } from './modules/anndata_overview'
include { SPECIFICITY } from './modules/specificity'
include { PLOT_SPECIFICITY } from './modules/plot_specificity'
include { SET_ANTIGENS_TYPE } from './modules/set_antigens_type'
include { PLOT_ANTIGENS_TYPE } from './modules/plot_antigens_type'


workflow {

    input_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )

    ANNDATA_OVERVIEW(input_ch)

    if (params.specificity == 'hpa') {

        SPECIFICITY(input_ch)

        PLOT_SPECIFICITY(
            SPECIFICITY.out.summary
        )

    } else {
        error "Unsupported specificity method: ${params.specificity}"
    }

    SET_ANTIGENS_TYPE(
        SPECIFICITY.out.classified
    )

    PLOT_ANTIGENS_TYPE(
        SET_ANTIGENS_TYPE.out.classified_labeled
    )
}
