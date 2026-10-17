/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { DECONVOLVE } from '../subworkflows/local/deconvolve'

include { ANALYZE } from '../subworkflows/local/analyze'

include  { QC } from '../modules/local/util/'

include { LOAD_DATASET } from '../subworkflows/local/load_dataset'

include { STAPLE_XSAMPLE } from '../modules/local/xsample/'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MULTIQC } from '../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow STAPLE {
    versions = channel.empty()

    // Load input paths and metadata
    INPUT_CHECK (
        file(params.input)
    )
    ch_multiqc_files = INPUT_CHECK.out.samplesheet_valid

    // Grab datasets
    ch_datasets = INPUT_CHECK.out.datasets
    LOAD_DATASET(ch_datasets)
    ch_adata = LOAD_DATASET.out.ch_adata
    ch_scrna = LOAD_DATASET.out.ch_scrna
    ch_matched_adata = LOAD_DATASET.out.ch_matched_adata
    versions = versions.mix(LOAD_DATASET.out.versions)

    // Report read data to MultiQC
    ch_report = LOAD_DATASET.out.ch_report

    // Deconvolve / cell type / use external
    DECONVOLVE (ch_datasets, ch_matched_adata, ch_scrna, ch_adata)
    versions = versions.mix(DECONVOLVE.out.versions)

    // Report deconvolution results to MultiQC
    ch_report = ch_report.mix( DECONVOLVE.out.ch_deconvolved.map {it -> [it.meta, it.obj, 'adata_counts']} ) //cell type counts
    ch_report = ch_report.mix( DECONVOLVE.out.ch_deconvolved.map {it -> [it.meta, it.obj, 'cell_probs']} )   //cell probs report

    // Filter deconvolved results to get anndata objects for squidpy spatial plots and ligand-receptor analysis
    ch_squidpy = DECONVOLVE.out.ch_deconvolved
                        .filter { it -> it.obj != null }
                        .filter { it -> it.obj.name.endsWith('.h5ad') }
                        .map { it-> [it.meta, it.obj] }

    // Analyze
    ANALYZE ( ch_squidpy )
    versions = versions.mix(ANALYZE.out.versions)

    // Cross-sample analysis
    if ( params.analyze.xsample ) {
        xsample_inputs = ANALYZE.out.adata.map{ it -> it[1] }.collect()
        STAPLE_XSAMPLE (xsample_inputs)
        versions = versions.mix(STAPLE_XSAMPLE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(STAPLE_XSAMPLE.out.multiqc_files)
    }

    // make reports
    QC( ch_report.filter { it -> it[1] != null && it[1] != '' && it[1] != [] }
            .map { it -> tuple(it[0], it[1], it[2]) } )
    ch_multiqc_files = ch_multiqc_files.mix(QC.out.report)
    versions = versions.mix(QC.out.versions)

    //collate versions
    softwareVersionsToYAML(versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'all_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set{ all_versions }
    ch_multiqc_files = ch_multiqc_files.mix(all_versions)

    // MultiQC
    // NOTE - will fail to find spaceranger reports unless the full path is provided
    // multiqc does not find spaceranger report for VisiumHD, address with #24
    def multiqc_config = file(params.multiqc_config).exists() ? file(params.multiqc_config) : []
    MULTIQC(
        ch_multiqc_files
            .collect()
            .ifEmpty([])
            .map { files -> [[id: 'multiqc'], files, multiqc_config, [], [], []] }
    )
    multiqc_report = MULTIQC.out.report.map { it -> it[1] }.toList()


    //emit:
      multiqc_report = multiqc_report           // channel: /path/to/multiqc_report.html
      versions       = versions                 // channel: [ versions.yml ]


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
