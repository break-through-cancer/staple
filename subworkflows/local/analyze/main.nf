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
        // if squidpy not run, just pass input adata along
        STAPLE_ATTACH_LIGREC(SQUIDPY_SPATIAL.out.adata.join(ligrec))
        attach_imscores_to = STAPLE_ATTACH_LIGREC.out.adata
        adata = attach_imscores_to
    } else {
        attach_imscores_to = SQUIDPY_SPATIAL.out.adata
        adata = attach_imscores_to
    }


    if (params.analyze.spacemarkers){
        STAPLE_ATTACH_IMSCORES(attach_imscores_to.join(imscores))
        attach_lrscores_to = STAPLE_ATTACH_IMSCORES.out.adata
        adata = attach_lrscores_to
        if (params.visium_hd){
            // in Visium HD, some spots have no reads, so make sure to keep all spots
            STAPLE_ATTACH_LRSCORES(attach_lrscores_to.join(lrscores))
            adata = STAPLE_ATTACH_LRSCORES.out.adata
        }
    }

    emit:
        versions
        adata
}