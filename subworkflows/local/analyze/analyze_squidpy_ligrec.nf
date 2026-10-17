include { SQUIDPY_LIGREC_ANALYSIS } from '../../../modules/local/squidpy/main'

workflow SQUIDPY_LIGREC {
    take:
        ch_squidpy
        ch_omnipath

    main:
        versions = channel.empty()
        SQUIDPY_LIGREC_ANALYSIS( ch_squidpy, ch_omnipath )
        versions = versions.mix(SQUIDPY_LIGREC_ANALYSIS.out.versions)

        ligrec = SQUIDPY_LIGREC_ANALYSIS.out.ligrec_interactions

    emit:
        versions
        ligrec
}