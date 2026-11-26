include { HARMONIZE_CELL_TYPES } from '../../../modules/local/util/'
include { STAPLE_ATTACH_HARMONIZE } from '../../../modules/local/util/'

workflow HARMONIZE {
    take:
        ch_deconvolved  // [meta, cell_probs, obj] for CoGAPS
        ch_patterns     // [meta, feature-by-pattern CSV]
    
    main:
    
    versions = channel.empty()
    
    ch_pattern_batch = ch_patterns
        .map { meta, patterns -> [meta.id.toString(), patterns] }
        .toList()
        .map { pairs -> tuple(pairs.collect { it[0] }, pairs.collect { it[1] }) }

    log.info "🔗 Harmonizing cell types across samples..."
    
    HARMONIZE_CELL_TYPES(
        ch_pattern_batch,
        params.harmonize.cut_height,
        params.harmonize.min_samples_per_cluster,
        params.harmonize.confidence_threshold
    )
    versions = versions.mix(HARMONIZE_CELL_TYPES.out.versions)
    
    // Map harmonized labels back to individual adata objects
    // For each sample's adata, join with the harmonization mapping file
    ch_with_mapping = ch_deconvolved
        .combine(HARMONIZE_CELL_TYPES.out.mapping)
        .map { item, mapping ->
            tuple(item.meta, item.obj, mapping, "cogaps")
        }
    
    STAPLE_ATTACH_HARMONIZE(ch_with_mapping)
    versions = versions.mix(STAPLE_ATTACH_HARMONIZE.out.versions)
    
    // Emit updated deconvolved channel with harmonized adatas
    // Format: [meta, cell_probs, obj] (same as input)
    ch_harmonized = ch_deconvolved
        .map { item -> tuple(item.meta, item.cell_probs) }
        .join(STAPLE_ATTACH_HARMONIZE.out.adata)
        .map { meta, cell_probs, adata -> [meta: meta, cell_probs: cell_probs, obj: adata] }
    
    emit:
        ch_harmonized
        mapping = HARMONIZE_CELL_TYPES.out.mapping
        consensus = HARMONIZE_CELL_TYPES.out.consensus
        diagnostics = HARMONIZE_CELL_TYPES.out.diagnostics
        dendrogram = HARMONIZE_CELL_TYPES.out.dendrogram
        versions
}
