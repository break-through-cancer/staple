# STAPLE: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

### Hardware requirements
Enough RAM and cpu to run spatial transcriptomics datasets that you have.

### Software requirements
Linux/mac OS with Docker, a slurm cluster with Docker or Singularity/Apptainer, or Cirro (https://cirro.bio/) is required to run STAPLE.

### STAPLE has been tested on
```
macOS 26.3.1 (a)
Cirro web platform (https://cirro.bio/)
slurm 23.11.6
github Ubuntu runners
```

## Just trying
 Run with test data in under 10 minutes on your laptop (needs docker and nextflow installed)! Use [this link](https://download-directory.github.io/?url=https://github.com/break-through-cancer/staple/tree/main/tests) to download ~30Mb of test data, then use the below commands in the terminal. Note that running STAPLE for the first time will trigger the download of the docker images, thus the actual time to run will be longer and will depend on the speed of your Internet connection. It may get up to 1h on slower speeds.

```
# make a clean directory called tests
mkdir tests

# extract the downloaded archive
unzip ~/Downloads/break-through-cancer\ staple\ main\ tests.zip -d tests

# navigate to tests/data
cd tests/data/samplesheets

# check that --max_memory and --max_cpus match your resources, run
nextflow run https://github.com/break-through-cancer/staple \
  --input multisample-test.csv \
  --max_memory 8GB \
  --max_cpus 4 \
  --outdir outs \
  -profile docker
```

Example terminal output:
```
 N E X T F L O W   ~  version 25.10.2

Launching `https://github.com/break-through-cancer/staple` [silly_rosalind] DSL2 - revision: fbf88bd923 [main]

executor >  local (28)
[7e/6aa75d] BTC:STAPLE:INPUT_CHECK:SAMPLESHEET_CHECK (multisample-test.csv) [100%] 1 of 1 ✔
[fb/4e7ece] BTC:STAPLE:LOAD_DATASET:ADATA_FROM_VISIUM (1)                   [100%] 2 of 2 ✔
[e6/388ec7] BTC:STAPLE:LOAD_DATASET:ADATA_ADD_METADATA (sample1)            [100%] 2 of 2 ✔
[11/078f2b] BTC:STAPLE:LOAD_DATASET:ATLAS_MATCH (sample1)                   [100%] 2 of 2 ✔
[70/e190b7] BTC:STAPLE:DECONVOLVE:RCTD (sample1)                            [100%] 2 of 2 ✔
[48/59abfc] BTC:STAPLE:DECONVOLVE:RCTD_PROBS (sample2)                      [100%] 2 of 2 ✔
[8b/23b9de] BTC:STAPLE:ANALYZE:SQUIDPY:SQUIDPY_LIGREC_ANALYSIS (sample1)    [100%] 2 of 2 ✔
[db/00ca6f] BTC:STAPLE:ANALYZE:SQUIDPY_SPATIAL_PLOTS (sample1)              [100%] 2 of 2 ✔
[28/f961d1] BTC:STAPLE:ANALYZE:STAPLE_ATTACH_LIGREC (sample1)               [100%] 2 of 2 ✔
[57/3f8f75] BTC:STAPLE:QC (9)                                               [100%] 10 of 10 ✔
[76/d2f213] BTC:STAPLE:MULTIQC                                              [100%] 1 of 1 ✔
Completed at: 01-Feb-2026 12:32:04
Duration    : 6m 54s
CPU hours   : 0.9
Succeeded   : 28
```
Examine the outputs in the `outs/` folder:
```
drwxr-xr-x@ 4 user  staff  128 Feb  1 12:25 adata
drwxr-xr-x@ 4 user  staff  128 Feb  1 12:26 atlas
drwxr-xr-x@ 4 user  staff  128 Feb  1 12:32 multiqc
drwxr-xr-x@ 8 user  staff  256 Feb  1 12:32 pipeline_info
drwxr-xr-x@ 4 user  staff  128 Feb  1 12:30 rctd
drwxr-xr-x@ 4 user  staff  128 Feb  1 12:31 squidpy
drwxr-xr-x@ 4 user  staff  128 Feb  1 12:32 staple
```

## Regular use
First, prepare a samplesheet with your input data that looks as follows: [samplesheet.csv](../samplesheet.csv), where each row represents a spatial transcriptomics sample.

The default named columns are following:

  * `sample` required, a unique identifier for the sample

  * `data_directory` required, path to the 10x spaceranger `outs` directory

  * `expression_profile`: optional, (leave blank if not using), reference expression profiles from matched scRNA (different local path to each scRNA atlas) or a local scRNA atlas (same path for each sample). If using a remotely stored atlas (such as CellXGene), rather pass `params.ref_scrna` and the atlas will be downloaded from the web.

Any extra columns will be treated as metadata and copied into the `meta` map, the resulting `.h5ad` object, and the final MultiQC report.


>[!IMPORTANT]
`export NXF_SINGULARITY_HOME_MOUNT=true` in order to allow matplotlib to write its logs (and avoid related error) if using singularity.

Run on Visium HD with RCTD (default) for cell typing and squidpy ligand-receptor analysis (default) using remote atlas annotation. In case of CellXGene atlas, the cell type column is always `cell_type`, so it does not need to be explicitly specified. In case a local `.h5ad` atlas is desired, specify the full path to it.
```bash
nextflow run break-through-cancer/staple \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --visium_hd <cell_segmentations/square_008um/square_016um/...> \
   --ref_scrna https://datasets.cellxgene.cziscience.com/d1d90d18-2109-412f-8dc0-e014e8abb338.h5ad
```
Run on Visium SD or HD with matched reference specified in the samplesheet and a custom cell type column `cell_type_column_name` in the reference:
```bash
nextflow run break-through-cancer/staple \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --ref_scrna_type_col cell_type_column_name
```

Pick a non-default reference-free deconvolution and ligand-receptor interaction tools:

```bash
nextflow run break-through-cancer/staple \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --deconvolve.<bayestme/cogaps/rctd> \
   --analyze.<spacemarkers/squidpy> \
```

Control squidpy [ligand-receptor interaction](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_ligrec.html) parameters:

Run in mouse data (or any other organism via ortholog conversion):
```bash
--sq_gr_ligrec_interactions_params '{"organism":"10090"}'
```

Use only CellPhoneDB:
```bash
--sq_gr_ligrec_interactions_params '{"resources": "CellPhoneDB"}'
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull break-through-cancer/staple
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [releases page](https://github.com/break-through-cancer/staple/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Arguments

### Pipeline

There is a number of arguments that may be passed to the pipeline to control its behavior and configuration. The list of all arguments is available in the `nextflow.config` or the nextflow_schema.json files, below we list the ones that we think are most relevant for typical usage.

#### Input format selection
Use `--visium`, `--visium_hd <table>`, `--xenium`, or `--anndata` depending on the format of your input data. For `anndata`, provide paths to `.h5ad` files; otherwise provide paths to the appropriate raw data directories. For Visium HD, select the associated table with, for example, `--visium_hd square_008um`.

#### Deconvolution / cell typing
RCTD drops cells/bins with less than 100 UMIs. Control that with `--umi_min <number>`. Setting this parameter to a lower value will include more cells/bins with fewer UMIs, while setting it to a higher value will be more stringent.

Use `--deconvolve.n_top_genes <number>` to control the number of top genes considered for deconvolution/cell typing. The default is 99999; setting it to a lower value can speed up computation and reduce memory usage.

#### Squidpy
Squidpy provides a lot of controls over database options for ligand-receptor analysis. Use `sq_gr_ligrec_interactions_params` to pass a JSON string with the desired parameters, for example: ` '{"resources": "CellPhoneDB"}'` to select only CellPhoneDB as the resource, or `'{"organism":"10090"}'` to select mouse instead of human.

The default way to estimate cell type co-occurrences is to split the whole slide into 50 bins. Change this, for example by using defined intervals like this: `sq_gr_co_occurrence_interval = '100,200,300'` or any other comma-separated list of bin sizes, or split in e.g. 10 bins instead: `sq_gr_co_occurrence_interval = 10`.

#### Cross-sample analysis
STAPLE performs cross-sample contrasts for metadata variables that are categorical and have exactly two observed levels across samples. Make sure to include relevant metadata columns in your samplesheet to enable meaningful cross-sample analysis. When performing such analysis, multiple testing correction is applied to control for false discoveries, but given the combination of many tests and potentially small sample sizes, the effects can get masked. Use `--analyze.spotlight` to specify which cell types to focus on in the cross-sample analysis, this reduces the number of tests and can help highlight significant effects. Example `--analyze.spotlight tumor`. Naturally, `'tumor'` should be a valid cell type present in your data.

Control how many top features are shown in the cross-sample analysis using the `--analyze.show_top` parameter. The default is 50. This affects only display and does not change the underlying data files produced and stored in `staple/reports/`.

`--analyze.only_spatial` flag will restrict the cross-sample analysis to only spatially variable genes, ignoring non-spatially variable genes in the pseudobulk analysis and ligand-receptor analysis.


### Core

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `...`

We recommend creating a profile specific to your own compute environment to ensure optimal performance and reproducibility, especially if running on a high-performance computing cluster or other specialized infrastructure.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.


## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
