
process SQUIDPY_SPATIAL {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:47a5a7292df74c7d4446609a5fae9676235292a60a1fa46ff02762cb3d10d0dc'

    input:
        tuple val(meta), path(adata)
    output:
        tuple val(meta), path("${prefix}/figures/spatial_scatter.png"),     emit: spatial_scatter_plot
        tuple val(meta), path("${prefix}/figures/interaction_matrix.png"),  emit: interaction_matrix_plot
        tuple val(meta), path("${prefix}/figures/co_occurrence*.png"),      emit: co_occurrence_plot
        tuple val(meta), path("${prefix}/figures/nhood_enrichment.png"),    emit: nhood_enrichment_plot
        tuple val(meta), path("${prefix}/figures/centrality_scores.png"),   emit: centrality_scores_plot
        tuple val(meta), path("${prefix}/squidpy.h5ad"),                    emit: adata
        path "versions.yml",                                                emit: versions


    script:
    // we may want to run for multiple deconvolution results and same sample
    source = adata.simpleName
    sample = "${meta.id}"
    prefix = task.ext.prefix ?: "${sample}/${source}"
    template 'spatial.py'
}

process SQUIDPY_LIGREC_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/break-through-cancer/btc-containers/scverse@sha256:47a5a7292df74c7d4446609a5fae9676235292a60a1fa46ff02762cb3d10d0dc'

    input:
        tuple val(meta), path(adata)
    output:
        tuple val(meta), path("${prefix}/**.png"),                            emit: ligrec_plots, optional: true
        tuple val(meta), path("${prefix}/**ligrec_interactions.pickle"),      emit: ligrec_interactions, optional: true
        tuple val(meta), path("${prefix}/**ligrec_metadata.csv"),             emit: ligrec_metadata, optional: true
        path "versions.yml",                                                  emit: versions

    script:
    // we may want to run for multiple deconvolution results and same sample
    source = adata.simpleName
    sample = "${meta.id}"
    prefix = task.ext.prefix ?: "${sample}/${source}"
    template 'ligrec.py'
}
