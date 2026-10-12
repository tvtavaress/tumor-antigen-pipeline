nextflow.enable.dsl=2

include { ANNDATA_OVERVIEW } from './modules/anndata_overview'
include { SPECIFICITY } from './modules/specificity'
include { PLOT_SPECIFICITY } from './modules/plot_specificity'
include { SET_ANTIGENS_TYPE } from './modules/set_antigens_type'
include { PLOT_ANTIGENS_TYPE } from './modules/plot_antigens_type'
include { NON_MALIGNANTS_HEALTHY_TISSUES } from './modules/non_malignants_healthy_tissues'
include { CANDIDATES_PER_TUMORS } from './modules/candidates_per_tumors'
include { PLOT_CANDIDATES_PER_TUMORS } from './modules/plot_candidates_per_tumors'


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

    NON_MALIGNANTS_HEALTHY_TISSUES(
        SET_ANTIGENS_TYPE.out.classified_labeled
    )

    if (params.pan_cancer) {

        gene_name_mapping_ch = params.gene_name_mapping
            ? Channel.value(file(params.gene_name_mapping, checkIfExists: true))
            : Channel.value(null)

        CANDIDATES_PER_TUMORS(
            input_ch,
            gene_name_mapping_ch,
            SPECIFICITY.out.classified
        )

        PLOT_CANDIDATES_PER_TUMORS(
            CANDIDATES_PER_TUMORS.out.binary_matrix
        )
    }
}