include { SPACEMARKERS } from './analyze_spacemarkers.nf'
include { SQUIDPY_LIGREC } from './analyze_squidpy_ligrec.nf'
include { SQUIDPY_SPATIAL } from '../../../modules/local/squidpy/main'
include { STAPLE_ATTACH_LIGREC } from '../../../modules/local/util/main'
include { STAPLE_ATTACH_LIGREC as STAPLE_ATTACH_IMSCORES } from '../../../modules/local/util/main'
include { STAPLE_ATTACH_LIGREC as STAPLE_ATTACH_LRSCORES } from '../../../modules/local/util/main'
include { OMNIPATH_GET } from '../../../modules/local/util/main'


workflow ANALYZE {

    take:
        ch_squidpy   // [meta, h5ad] from DECONVOLVE
    main:

    versions = channel.empty()
    ligrec = channel.empty()
    imscores = channel.empty()
    lrscores = channel.empty()
    ligrec_db = channel.empty()

    // get interactions if needed: param suplied or get from omnipath
    if (params.analyze.spacemarkers || params.analyze.squidpy){
        if (params.ref_ligrec) {
            ligrec_db = file(params.ref_ligrec)
        } else {
            OMNIPATH_GET()
            versions = versions.mix(OMNIPATH_GET.out.versions)
            ligrec_db = OMNIPATH_GET.out.db
        }
    }

    // ligrec - spacemarkers if requested
    if (params.analyze.spacemarkers){
        SPACEMARKERS(ch_squidpy, ligrec_db)
        versions = versions.mix(SPACEMARKERS.out.versions)
        // pass ligrec results along
        imscores = imscores.mix(SPACEMARKERS.out.imscores)
        lrscores = lrscores.mix(SPACEMARKERS.out.lrscores)
    }

    // ligrec - squidpy if requested
    if (params.analyze.squidpy){
        SQUIDPY_LIGREC( ch_squidpy, ligrec_db )
        versions = versions.mix(SQUIDPY_LIGREC.out.versions)
        ligrec = ligrec.mix(SQUIDPY_LIGREC.out.ligrec)
    }

    // do basic analysis anyway
    SQUIDPY_SPATIAL( ch_squidpy )
    versions = versions.mix(SQUIDPY_SPATIAL.out.versions)

    // wrap up - collect results from tools and save
    // TODO: rewrite to collect ligrecs and join once
    if (params.analyze.squidpy){
        // remainder keeps samples without ligrec; [] stages no file so staple.h5ad is still written
        ch_with_ligrec = SQUIDPY_SPATIAL.out.adata.join(ligrec, remainder: true)
            .map { meta, ad, lr -> [meta, ad, lr ?: []] }
        STAPLE_ATTACH_LIGREC(ch_with_ligrec)
        attach_imscores_to = STAPLE_ATTACH_LIGREC.out.adata
        adata = attach_imscores_to
    } else {
        attach_imscores_to = SQUIDPY_SPATIAL.out.adata
        adata = attach_imscores_to
    }


    if (params.analyze.spacemarkers){
        ch_with_im = attach_imscores_to.join(imscores, remainder: true)
            .map { meta, ad, im -> [meta, ad, im ?: []] }
        STAPLE_ATTACH_IMSCORES(ch_with_im)
        attach_lrscores_to = STAPLE_ATTACH_IMSCORES.out.adata
        ch_with_lr = attach_lrscores_to.join(lrscores, remainder: true)
            .map { meta, ad, lr -> [meta, ad, lr ?: []] }
        STAPLE_ATTACH_LRSCORES(ch_with_lr)
        adata = STAPLE_ATTACH_LRSCORES.out.adata
    }

    emit:
        versions
        adata
}