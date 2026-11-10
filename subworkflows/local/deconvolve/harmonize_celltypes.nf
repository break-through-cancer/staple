include { HARMONIZE_CELL_TYPES } from '../../../modules/local/util/'
include { ATTACH_HARMONIZED_LABELS } from '../../../modules/local/util/'

workflow HARMONIZE {
    take:
        ch_deconvolved  // [meta, cell_probs, obj_adata] from DECONVOLVE
    
    main:
    
    versions = channel.empty()
    
    // Extract reference-free deconvolution only (cogaps or bayestme with cell_probs)
    // These have cell_probs CSV output; RCTD and external don't go through harmonization
    ch_ref_free = ch_deconvolved
        .filter { it.cell_probs != null }
        .map { it -> [it.meta, it.cell_probs, it.obj] }
    
    // Collect as strings to avoid Nextflow path collision detection
    // Convert path objects to their string paths
    ch_file_paths_str = ch_ref_free
        .map { meta, cell_probs, adata -> cell_probs.toString() }
        .collect()
        .map { list -> list.join('\n') }  // Join with newlines for file list
    
    ch_sample_ids_list = ch_ref_free
        .map { meta, cell_probs, adata -> meta.id }
        .collect()
        .map { it.join(',') }
    
    log.info "🔗 Harmonizing cell types across samples..."
    
    // Call harmonization with file paths as a string
    HARMONIZE_CELL_TYPES(
        ch_file_paths_str,
        ch_sample_ids_list,
        params.harmonize.cut_height,
        params.harmonize.min_samples_per_cluster,
        params.harmonize.confidence_threshold
    )
    versions = versions.mix(HARMONIZE_CELL_TYPES.out.versions)
    
    // Map harmonized labels back to individual adata objects
    // For each sample's adata, join with the harmonization mapping file
    ch_with_mapping = ch_ref_free
        .combine(HARMONIZE_CELL_TYPES.out.mapping)
        .map { meta, cell_probs, adata, mapping -> 
            tuple(meta, adata, mapping, "cogaps")
        }
    
    ATTACH_HARMONIZED_LABELS(ch_with_mapping)
    versions = versions.mix(ATTACH_HARMONIZED_LABELS.out.versions)
    
    // Emit updated deconvolved channel with harmonized adatas
    // Format: [meta, cell_probs, obj] (same as input)
    ch_harmonized = ATTACH_HARMONIZED_LABELS.out.adata
        .map { meta, adata -> 
            [meta: meta, cell_probs: null, obj: adata]
        }
    
    emit:
        ch_harmonized
        mapping = HARMONIZE_CELL_TYPES.out.mapping
        consensus = HARMONIZE_CELL_TYPES.out.consensus
        diagnostics = HARMONIZE_CELL_TYPES.out.diagnostics
        dendrogram = HARMONIZE_CELL_TYPES.out.dendrogram
        versions
}
