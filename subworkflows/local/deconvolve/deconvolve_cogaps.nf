
include { COGAPS_ADATA2DGC;
          COGAPS as COGAPS_MAIN;
          COGAPS_PREPROCESS } from '../../../modules/local/cogaps/main'
include { CELL_TYPES_FROM_COGAPS as COGAPS_CELL_TYPES } from '../../../modules/local/util/'

workflow COGAPS {
    take:
        ch_adata // channel of [meta, adata]

    main:

    def patterns = params.npatterns.split(',').collect { it -> it.toInteger() }
    ch_patterns = channel.from(patterns)
    versions = channel.empty()

    //example channel with cparams
    ch_fixed_params = channel.of([seed: params.seed, niterations: params.niterations, sparse: params.sparse, distributed: params.distributed, nsets:params.nsets, nthreads:1, asynchronous_updates:params.asynchronous_updates])

    ch_cparams = ch_patterns
      .combine(ch_fixed_params)
      .map { it -> tuple([id:it[0].toString(), npatterns:it[0], seed:it[1].seed, niterations:it[1].niterations, sparse:it[1].sparse, distributed:it[1].distributed, nsets:it[1].nsets, nthreads:it[1].nthreads, asynchronous_updates:it[1].asynchronous_updates]) }

    // convert adata to dgCMatrix
    COGAPS_ADATA2DGC(ch_adata)
    versions = versions.mix(COGAPS_ADATA2DGC.out.versions)


    // preprocess dgCMatrix
    ch_preprocess = COGAPS_ADATA2DGC.out.dgCMatrix
      .map { it -> tuple(it[0], it[1]) }

    COGAPS_PREPROCESS(ch_preprocess)
    versions = versions.mix(COGAPS_PREPROCESS.out.versions)

    ch_input = COGAPS_PREPROCESS.out.dgCMatrix
      .map { it -> tuple(it[0], it[1]) }
      .combine(ch_cparams)

    COGAPS_MAIN(ch_input)
    versions = versions.mix(COGAPS_MAIN.out.versions)
    ch_cogaps = COGAPS_MAIN.out.cogapsResult.map { it -> tuple(it[0], it[1]) }

    COGAPS_CELL_TYPES(ch_cogaps)
    versions = versions.mix(COGAPS_CELL_TYPES.out.versions)
    ch_deconvolved = COGAPS_CELL_TYPES.out.cogaps_cell_types.map { it -> tuple(it[0], it[1]) }

    emit:
        ch_deconvolved
        versions
}