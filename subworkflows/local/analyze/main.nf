include { SPACEMARKERS } from './analyze_spacemarkers.nf'
include { SQUIDPY_LIGREC } from './analyze_squidpy_ligrec.nf'
include { SQUIDPY_SPATIAL } from '../../../modules/local/squidpy/main'
include { STAPLE_ATTACH_LIGREC } from '../../../modules/local/util/main'
include { STAPLE_ATTACH_LIGREC as STAPLE_ATTACH_IMSCORES } from '../../../modules/local/util/main'
include { STAPLE_ATTACH_LIGREC as STAPLE_ATTACH_LRSCORES } from '../../../modules/local/util/main'


workflow ANALYZE {

    take: 
        ch_sm_inputs   // from DECONVOLVE
        ch_squidpy
    main:

    versions = channel.empty()
    ligrec = channel.empty()
    imscores = channel.empty()
    lrscores = channel.empty()


    // ligrec - spacemarkers if requested
    if (params.analyze.spacemarkers){
        SPACEMARKERS(ch_sm_inputs)
        versions = versions.mix(SPACEMARKERS.out.versions)
        // pass ligrec results along
        imscores = imscores.mix(SPACEMARKERS.out.imscores)
        lrscores = lrscores.mix(SPACEMARKERS.out.lrscores)
    }

    // ligrec - squidpy if requested
    if (params.analyze.squidpy){
        SQUIDPY_LIGREC( ch_squidpy )
        versions = versions.mix(SQUIDPY_LIGREC.out.versions)
        ligrec = ligrec.mix(SQUIDPY_LIGREC.out.ligrec)
    }

    // do basic analysis anyway
    SQUIDPY_SPATIAL( ch_squidpy )
    versions = versions.mix(SQUIDPY_SPATIAL.out.versions)

    // wrap up - collect results from tools and save
    // TODO: rewrite to collect ligrecs and join once
    if (params.analyze.squidpy){
        // remainder: keep samples whose ligrec output is missing (optional output / failed ignored)
        ch_with_ligrec = SQUIDPY_SPATIAL.out.adata.join(ligrec, remainder: true)
        STAPLE_ATTACH_LIGREC(ch_with_ligrec.filter { it[2] != null })
        // samples without ligrec pass through unchanged
        attach_imscores_to = STAPLE_ATTACH_LIGREC.out.adata.mix(
            ch_with_ligrec.filter { it-> it[2] == null }.map { meta, ad, _lr -> [meta, ad] }
        )
        adata = attach_imscores_to
    } else {
        attach_imscores_to = SQUIDPY_SPATIAL.out.adata
        adata = attach_imscores_to
    }


    if (params.analyze.spacemarkers){
        ch_with_im = attach_imscores_to.join(imscores, remainder: true)
        STAPLE_ATTACH_IMSCORES(ch_with_im.filter { it[2] != null })
        attach_lrscores_to = STAPLE_ATTACH_IMSCORES.out.adata.mix(
            ch_with_im.filter { it -> it[2] == null }.map { meta, ad, _im -> [meta, ad] }
        )
        adata = attach_lrscores_to
        if (params.visium_hd){
            ch_with_lr = attach_lrscores_to.join(lrscores, remainder: true)
            STAPLE_ATTACH_LRSCORES(ch_with_lr.filter { it -> it[2] != null })
            adata = STAPLE_ATTACH_LRSCORES.out.adata.mix(
                ch_with_lr.filter { it -> it[2] == null }.map { meta, ad, _lr -> [meta, ad] }
            )
        }
    }

    emit:
        versions
        adata
}