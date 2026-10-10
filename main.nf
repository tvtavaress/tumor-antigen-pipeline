nextflow.enable.dsl=2

include { ANNDATA_OVERVIEW } from './modules/anndata_overview'

workflow {
    input_ch = Channel.fromPath(
        params.input,
        checkIfExists: true
    )

    ANNDATA_OVERVIEW(input_ch)
}
