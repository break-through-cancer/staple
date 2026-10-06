
include { BAYESTME } from './deconvolve_bayestme'
include { COGAPS } from './deconvolve_cogaps'
include { RCTD } from '../../../modules/local/rctd/rctd'
include { ATTACH_CELL_PROBS as EXTERNAL_PROBS } from '../../../modules/local/util/'
include { ATTACH_CELL_PROBS as RCTD_PROBS } from '../../../modules/local/util/'
include { ATTACH_CELL_PROBS as COGAPS_PROBS } from '../../../modules/local/util/'

workflow DECONVOLVE {
    take:
        ch_datasets          // from LOAD_DATASET, for BayesTME, as it has it's own Visium reader
        ch_matched_adata     // from LOAD_DATASET, for RCTD and any reference based deconvolution
        ch_scrna             // from LOAD_DATASET, for RCTD and any reference based deconvolution
        ch_adata             // from LOAD_DATASET, for CoGAPS and any reference free deconvolution
    main:

    versions = channel.empty()
    ch_deconvolved = channel.empty()  //outputs: meta, optional cell type probs, optional deconvolution object

    // Skip deconvolution and work with existing cell types
    // cell_type column expected as per cellXgene specification
    if (params.deconvolve.skip) {
        ch_deconvolved = ch_deconvolved.mix(ch_adata.map { it -> [meta:it[0], cell_probs:null, obj:it[1]] })
    }

    // Grab external deconvolution results
    // CODA annotation channel - or any other external csv annotaion
    if (params.deconvolve.external && !(params.deconvolve.skip || params.deconvolve.skip == 'true')) {
        ch_external = ch_datasets.map { it -> tuple(it.meta, it.data_directory) }
            .flatMap { item -> 
                def meta = item[0]
                def data_path = item[1]
                def coda_files = []
                data_path.eachFileRecurse { file ->
                    if (file.name.endsWith('cellular_compositions.csv')) {
                        coda_files.add(file)
                    }
                }
                coda_files.collect { file -> [meta, file] }
            }
        //join external cell probs to adata
        ch_cell_probs_input = ch_external.join( ch_adata ).map { it -> tuple(it[0], it[1], it[2], "external") }
        EXTERNAL_PROBS(ch_cell_probs_input)
        ch_deconvolved = ch_deconvolved.mix(ch_external.join(EXTERNAL_PROBS.out.adata)
            .map { it -> [meta:it[0], cell_probs:it[1], obj:it[2]] })
    }

    // BayestME deconvolution and plots, run only if not hd as the tool does not support it
    if(!params.visium_hd && params.deconvolve.bayestme && !(params.deconvolve.skip || params.deconvolve.skip == 'true')) {
        BAYESTME(ch_datasets)
        ch_deconvolved = ch_deconvolved.mix(BAYESTME.out.ch_deconvolved.map { it -> [meta:it[0], cell_probs:null, obj:it[1]] })
        versions = versions.mix(BAYESTME.out.versions)
    }

    // RCTD reference-based deconvolution and plots
    // plots are temporary as there is the idea to plot
    // all deconvolution stats with the same process/workflow
    if(params.deconvolve.rctd && !(params.deconvolve.skip || params.deconvolve.skip == 'true')) {
        ch_rctd_input = ch_scrna.join(ch_matched_adata)
        RCTD( ch_rctd_input )
        versions = versions.mix(RCTD.out.versions)
        ch_rctd_output = RCTD.out.rctd_cell_types

        // join rctd cell types to adata
        ch_cell_probs_input = ch_rctd_output.join( ch_adata ).map { it -> tuple(it[0], it[1], it[2], "rctd") }
        RCTD_PROBS(ch_cell_probs_input)
        ch_rctd_output = ch_rctd_output.join(RCTD_PROBS.out.adata)
            .map { it -> [meta:it[0], cell_probs:it[1], obj:it[2]] }
        ch_deconvolved = ch_deconvolved.mix(ch_rctd_output)
        versions = versions.mix(RCTD.out.versions)
    }

    //CoGAPS reference-free spatially unaware deconvolution
    if(params.deconvolve.cogaps && !(params.deconvolve.skip || params.deconvolve.skip == 'true')) {
        COGAPS( ch_adata )
        versions = versions.mix(COGAPS.out.versions)
        ch_cogaps_output = COGAPS.out.ch_deconvolved

        //join cogaps cell types to adata
        COGAPS_PROBS( ch_cogaps_output
                        .join( ch_adata ).map { it-> tuple(it[0], it[1], it[2], "cogaps") } ) //meta, cell types, adata, label
        ch_cogaps_output = ch_cogaps_output.join(COGAPS_PROBS.out.adata)
            .map { it-> [meta:it[0], cell_probs:it[1], obj:it[2]] }
        ch_deconvolved = ch_deconvolved.mix(ch_cogaps_output)
        versions = versions.mix(COGAPS.out.versions)
    }

    emit:
        ch_deconvolved
        versions
}