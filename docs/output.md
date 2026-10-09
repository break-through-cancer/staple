# STAPLE: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Load dataset](#load-dataset) - Load spatial dataset and optionally the single-cell reference dataset.
- [Deconvolve](#deconvolve) - Perform deconvolution / cell type annotation of the spatial dataset.
- [Analyze](#analyze) - Perform downstream analyses on the spatial dataset, including spatial statistics and ligand-receptor analysis.
- [QC](#qc) - Perform quality control checks on the data.
- [Staple](#staple) - Perform cross-sample analysis.
- [MultiQC](#multiqc) - Generate a comprehensive report summarising all metrics.


### Load dataset

<details markdown="1">
<summary>Output files</summary>

- `adata/`
  - `<sample>/`
    - `adata.h5ad`: AnnData object containing the spatial dataset.
- `atlas/`
  - `<sample>/`
    - `adata.h5ad`: AnnData object matched to the single-cell reference dataset, if provided.
  - `atlas.h5ad`: AnnData object containing the single-cell reference dataset, if provided.

</details>

Load dataset uses [spatialdata.io](https://spatialdata.scverse.org/projects/io/en/latest/) to load the spatial dataset and optionally the single-cell reference dataset. The output is converted to AnnData objects and stored in the `adata/` directory. If a single-cell reference dataset is provided, the matched AnnData objects are stored in the `atlas/` directory and the atlas object is stored as `atlas.h5ad`.

### Deconvolve
<details markdown="1">
<summary>Output files</summary>

- `<tool>/`
  - `<sample>/`
    - `<tool>.h5ad`: AnnData object containing cell type annotation in `.obs`.

</details>

The deconvolution step performs cell type annotation of the spatial dataset using the provided single-cell reference dataset, by reference-free methods, or using an external annotation provided by the user. The output is stored in the `<tool>/` directory as an AnnData object containing the cell type annotation in `.obs`, where `<tool>` is the name of the deconvolution method used (e.g. `RCTD`, `CoGAPS`, `BayesTME`, or `external`). If a tool does not directly output an AnnData object, the tool-agnostic output will be saved, too.

### QC

QC module does not produce any output files, but the results of the QC checks are visualised in the MultiQC report and are also available in the `multiqc/multiqc_data/` directory.

### Analyze
<details markdown="1">
<summary>Output files</summary>

- `squidpy/`
  - `<sample>/`
      - `<tool>/`
        - `squidpy.h5ad`: object containing the spatial analyses in `.obs` and `.uns`.
        - `figures/`: directory containing spatial analysis static per-sample images
        - `ligrec/`
          - `figures/`: figures directory containing ligand-receptor analysis static per-sample images
          - `ligrec-interactions.pickle`: Pickle object containing the results of the ligand-receptor analysis.
- `spacemarkers/`
  - `<sample>/`
      - `<source>/`
        - `IMscores.rds`: SpaceMarkers interaction scores for genes (directed or undirected, see usage)
        - `LRscores.rds`: SpaceMarkers ligand-receptor interaction scores (directed mode only)

</details>

Spatial metrics are computed using [Squidpy](https://squidpy.readthedocs.io/en/stable/). Ligand-receptor analysis is done with Squidpy's ligand-receptor analysis functionality, and also with [SpaceMarkers](https://github.com/DeshpandeLab/SpaceMarkers). The tools provide different outputs that are saved as intermediate outputs in the `squidpy/` and `spacemarkers/` directories, respectively. Static images of the spatial analyses are also saved in the `figures/` directories. Downstream these results are integrated into the MultiQC report and the final AnnData objects.

### Staple
<details markdown="1">
<summary>Output files</summary>

- `staple/`
  - `<sample>/`
    - `staple.h5ad`: AnnData object containing all results in `.obs` and `.uns`.
  - `reports/`
    - `mqc/`
      - `ligrec_diff_response_results.json`: results of the ligand-receptor differential response analysis.
      - `neighbors_mqc.json`: results of the neighborhood analysis.
      - `dese2_diff_response_results.json`: results of the DESeq2 differential response analysis. All strata (cell types) are included in the same report, with the strata name included in the ID of each result.
    - `ligrec_diff_response_results.csv`: results of the ligand-receptor differential response analysis.
    - `Moran_I_diff_response_results.csv`: results of the Moran's I differential response analysis.
    - `deseq2_diff_response_results_<strata>.csv`: results of the DESeq2 differential response analysis for a specific strata (cell_type).
    - `pseudobulk.h5ad` - AnnData object containing the pseudobulked data used for the DESeq2 differential response analysis for downstream contrasts.

</details>
The `staple/` directory contains the final integrated results of the pipeline, including the final AnnData objects containing all results in `.obs` and `.uns`.


### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline results are visualised in the report and further statistics are available in the report data directory. MultiQC also provides means to chat about the results with LLMs via an API, or by exporting the prompt data from the report.

Results generated by MultiQC collate pipeline reports from supported tools. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
