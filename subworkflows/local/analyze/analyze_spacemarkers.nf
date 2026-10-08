include { SPACEMARKERS as SPACEMARKERS_RUN } from '../../../modules/local/spacemarkers/nextflow/main'
include { SPACEMARKERS_CHECK } from '../../../modules/local/util/'
include { SPACEMARKERS_HARMONIZE as SPACEMARKERS_HARMONIZE_IMSCORES } from '../../../modules/local/util/'
include { SPACEMARKERS_HARMONIZE as SPACEMARKERS_HARMONIZE_LRSCORES } from '../../../modules/local/util/'


workflow SPACEMARKERS {

    take:
        ch_adata   // [meta, h5ad] after cell typing
    main:
        versions = channel.empty()

        // keep only samples that have latent features in adata.uns
        SPACEMARKERS_CHECK( ch_adata )
        versions = versions.mix(SPACEMARKERS_CHECK.out.versions)
        ch_eligible = SPACEMARKERS_CHECK.out.checked
            .filter { _meta, _adata, eligible -> eligible.trim() == 'true' }
            .map { meta, adata, _eligible -> [meta, adata] }

        SPACEMARKERS_RUN( ch_eligible )
        // IMScores: IMscores.rds with row names holding gene names, one column per pattern pair
        // LRscores: LRscores.rds with ligand-receptor pairs as row names, directed mode only
        versions = versions.mix(SPACEMARKERS_RUN.out.versions)

    SPACEMARKERS_HARMONIZE_IMSCORES( SPACEMARKERS_RUN.out.IMscores )
    SPACEMARKERS_HARMONIZE_LRSCORES( SPACEMARKERS_RUN.out.LRscores )

    imscores = SPACEMARKERS_HARMONIZE_IMSCORES.out.spacemarkers
    lrscores = SPACEMARKERS_HARMONIZE_LRSCORES.out.spacemarkers

    emit:
        versions
        imscores
        lrscores

}