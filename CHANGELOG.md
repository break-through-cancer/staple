# btc/spatial: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.6.0 - 2026-10-26

### `Added`
- custom ligand-receptor interaction databases and MultiQC configuration
- OmniPath support for SpaceMarkers ligand-receptor analysis
- filter interactions db based on a query to any field

### `Fixed`
- SpaceMarkers input handling and empty ligand-receptor results
- t-test result formatting and MultiQC sample selection

## v2.5.0 - 2026-10-06

### `Added`
- option to skip deconvolution
- support for custom prompts and improved MultiQC reporting

## v2.4.2 - 2026-09-25

### `Added`
- improved custom prompting and cross-sample reporting

## v2.4.1 - 2026-09-14

### `Fixed`
- spotlight mode and cross-sample analysis issues

## v2.3.1 - 2026-08-28

### `Added`
- ignore self-neighbors by default in cross-sample analysis

## v2.3.0 - 2026-07-27

### `Added`
- support for local reference atlases and additional input formats
- ligand-receptor database selection and mitochondrial gene filtering
- expanded cross-sample analysis and MultiQC reporting

### `Fixed`
- RCTD doublet-mode fallback and neighbor reports for absent cell types

## v2.2.0 - 2025-11-23

### `Added`
- collect all outputs to final anndatas
- cross-sample analysis module
- extended multiqc reports
- extended tests 

## v2.1.0 - 2025-11-23

### `Fixed`
- moved all tool-related params samplesheet to params
- fixed all nextflow linter errors
- brought back MultiQC report
- harmonized all deconvolution method io

### `Added`
- segmentation support for Visium HD
- samplesheet with all metadata fields to MultiQC
- metamap extends with every column from samplesheet except data paths
- random seed to SpaceMarkers, RCTD, Squidpy modules
- NA as value option for cell_type, used downstream from deconvolution/cell typing
- basic QC for Anndata objects

### `Deprecated`
- removed samplesheet params `n_cell_types,n_top_genes,bleeding_correction,run_bayestme,run_cogaps,run_spacemarkers,find_annotations`
- disabled SpaceMarkers SD plots
- removed bleeding_correction Cirro param for simplicity


## v1.3.1 - 2025-09-06

### `Fixed`
- deprecated check_max that stopped working on cirro

## v1.3.0 - 2025-08-28

### `Added`
- squidpy ligand-receptor interactions (including plot of top 25)
- pipeline tests for visium SD and HD
- SpaceMarkers support for visium HD
- doublet_mode param for rctd

### `Fixed`
- collect versions from utils subworkflow (missing before)
- collect versions from all squidpy submodules
- RCTD module trying to use `adata.raw.X` when no `adata.X` present
- refactor config files
- refactor nextflow schema

## v1.2.0 - 2025-07-11

### `Added`
- optionally look for external (e.g. CODA) annotation in the input to run spaceMarkers on
- spatial overlap and interaction plots for spaceMarkers
- RCTD deconvolution with single cell atlas reference
- Visium HD support (RCTD)
- Squidpy module for plotting deconvolution results and spatial analysis

### `Fixed`
- refactor: separate subworkflow for data loading
- refactor: separate subworkflow for BayesTME to avoid it clashing with VisiumHD pipeline

### `Deprecated`
- spatial transcriptional programs from the BayesTME package
- BAYESTME_LOAD_SPACERANGER deprecated to allow Visium HD (and possible more) input
- do not reuse BayesTME preprocessing for CoGAPS


## v1.1.0 - 2024-12-16

### `Added`
- multiqc report for Spaceranger
- multiqc report for SpaceMarkers
- standalone SpaceMarkers module to run after any/both deconvolutions
- CoGAPS deconvolution module
- major analysis steps can now be turned on/off
- number of top genes to keep as parameter for deconvolution
- radically run CoGAPS after BayesTME preprocessing

### `Fixed`
- Cirro not being able to run on multiple samples
- Cirro showing 'sample' instead of the proper sample name
- SpaceMarkers returning no interactions in BayesTME results
- SpaceMarkers not running multithreaded
- spatial transcriptional programs running if enable = False


## v1.0dev - [date]

Initial release of btc/spatial, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
