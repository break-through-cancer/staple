#!/usr/bin/env python3
# Cross-sample analysis template script
# Input is a space-delimited string with adatas

import os
import logging
import pandas as pd
import scipy as sp
import anndata as ad
import scanpy as sc
import numpy as np
import json
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def ligrec_from_adatas(adatas, type='ligrec_means', axis=1,
                            samples=None, spotlight=None):

    # extract stat from each adata
    ligrecs = [x.uns[type] for x in adatas if type in x.uns]

    # combine sample level
    combined = pd.concat(ligrecs, keys=samples, axis=axis)

    # move cell types to columns and perform the t-test
    res = combined.stack(future_stack=True)

    # if a cell type or pair has been specified, filter to that only
    if spotlight:
        sp_index = [np.all([y in x for y in spotlight ]) for x in res.index.get_level_values(-1)]
        res = res[sp_index]

    return res

def heatmap_report(adatas, spotlight=None, groups=None, show=100, filter=0.05, tool=None):
    samples = [a.obs['id'].unique()[0] for a in adatas]
    pvalues = None
    if tool =='squidpy_ligrec':
        ligrecs = ligrec_from_adatas(adatas, type='ligrec_means', spotlight=spotlight, samples=samples)
        pvalues = ligrec_from_adatas(adatas, type='ligrec_pvalues', spotlight=spotlight, samples=samples)
    elif tool =='spacemarkers_LRscores':
        ligrecs = ligrec_from_adatas(adatas, type='LRscores', spotlight=spotlight, samples=samples)
    elif tool =='Moran_I':
        #Moran's I is not per cell type pair, so no spotlighting (yet)
        ligrecs = ligrec_from_adatas(adatas, type='moranI', spotlight=None, samples=samples)
        #pick only the Moran's I value index
        ligrecs = ligrecs[ligrecs.index.get_level_values(-1) == 'I']
        ligrecs = ligrecs.droplevel(-1)

    if pvalues is not None:
        # filter sample ligrec pairs by p-value to reduce multiple testing
        sig_mask = (pvalues <= filter).any(axis=1)
        ligrecs = ligrecs[sig_mask]

    #join multiindex of squidpy ligrec into single index
    if(isinstance(ligrecs.index, pd.MultiIndex)):
        ligrecs.index = ['-'.join(map(str, idx)) for idx in ligrecs.index]

    if groups is not None:
        ligrecs_ttest = xsample_ttest(ligrecs, groups[0], groups[1])
        ligrecs_ttest_sig = ligrecs_ttest[ligrecs_ttest['pval_adj'] <= filter]

        #mark significant with a star in the index for plotting purposes
        ligrecs_ttest.index = [idx + '*' if idx in ligrecs_ttest_sig.index \
            else idx for idx in ligrecs_ttest.index]

        if(len(ligrecs_ttest_sig) == 0):
            log.warning(f"No significant differential interactions found for {tool}.")
        else:
            log.info(f"{len(ligrecs_ttest_sig)} significant differential \
                interactions found for {tool} with p_adj<={filter}.")

        res = ligrecs_ttest.sort_values('pval', ascending=True)
        memo = f"Top {show} differential interactions across samples. There are \
                {len(ligrecs_ttest_sig['pval_adj'])} significant differential \
                interactions (p_adj<={filter}) found between \
                {groups[0]} and {groups[1]}, marked with *."

    else:
        ligrecs['mean'] = ligrecs.mean(axis=1)
        res = ligrecs.sort_values('mean', ascending=False)
        memo = f"Top {show} mean interactions across samples shown as no groups \
                 were specified in the sample sheet."

    res_show = res[samples][:show]
    res_dict = res_show.to_dict()

    mqc_report = {
        "id": f"{tool}_interactions",
        "description": memo,
        "plot_type": "heatmap",
        "pconfig": {
            "ylab": "Sample",
            "ycats_samples": True,
            "xcats_samples": False,
            "xlab": "Interaction",
            "zlab": "Score",
            "title": "Interaction scores",
            "square": False
        },
        "data": res_dict
    }
    return mqc_report, res


def pseudobulk_adatas(adatas, vars=None, only_spatial=False):
    # collect pseudobulk expression for each adata
    # if only_spatial, select only genes with spatially variable expression
    pb_adatas = []
    for adata in adatas:
        # check that there are no NA values in grouping
        if vars is not None:
            for var in vars:
                if adata.obs[var].isna().any():
                    log.warning(f"NA values found in grouping variable {var} for sample {adata.obs['id'].unique()[0]}.\
                        Replacing NA with 'NA' string for pseudobulk grouping.")
                    series = adata.obs[var]
                    if isinstance(series.dtype, pd.CategoricalDtype):
                        if 'NA' not in series.cat.categories:
                            series = series.cat.add_categories(['NA'])
                        adata.obs[var] = series.fillna('NA')
                    else:
                        adata.obs[var] = series.astype(object).fillna('NA')
        adata = adata.to_memory()
        if only_spatial:
            adata = adata[:, adata.var['spatially_variable']==True]
        pb_adatas.append(sc.get.aggregate(adata, by=vars, func='sum'))
    pb_adata = ad.concat(pb_adatas, axis=0, join='outer')
    counts = pb_adata.layers['sum']
    # keep raw counts in X for pydeseq2, handling both dense and sparse matrices
    if sp.sparse.issparse(counts):
        counts = counts.tocsr().copy()
        # replace NaNs with zeros in the underlying sparse data array
        data = counts.data
        nan_mask = np.isnan(data)
        if np.any(nan_mask):
            data[nan_mask] = 0
    else:
        counts = np.asarray(counts).copy()
        counts[np.isnan(counts)] = 0
    pb_adata.X = counts
    del pb_adata.layers['sum']

    return pb_adata


def pydeseq_results(pb_adata, spotlight=None, cpus=1, design=None, contrast=None):
    # estimate size factors and dispersions with pydeseq2
    inference = DefaultInference(n_cpus=cpus)
    dds = DeseqDataSet(adata=pb_adata,
                       design=design,
                       inference=inference)
    dds.deseq2()
    de = DeseqStats(dds, contrast=contrast)

    return de

def de_report(de_dict, spotlight=None, filter=0.05, show=100, contrast=None, p='padj'):
    #plot logfoldchangs vs -log padj for all genes
    # multiqc want this format: "gene1": {"x": 1, "y": 2}
    # x is log2 fold change, y is -log10 adjusted p-value
    # show_p is mainly for testng purposes to allow plotting
    reports = []
    for de in de_dict.values():
        res_sig = de[de[p] <= filter]
        res_sig_show = res_sig.sort_values(p).head(show)
        res_sig_show.rename(columns={"log2FoldChange": "x", p: "y"}, inplace=True)
        res_sig_show['y'] = -np.log10(res_sig_show['y'])
        res_dict = res_sig_show.loc[:, ['x','y']].to_dict(orient='index')
        reports.append(res_dict)
    memo = f"DESeq2 results for {contrast[0]} variable with significant DE genes ({p}<={filter}) between {contrast[1]} and {contrast[2]}. Log2 fold change (X) and -log10 adjusted p-value (Y) shown."
    mqc_report = {
        "id": f"deseq2_{contrast[0]}",
        "description": memo,
        "plot_type": "scatter",
        "pconfig": {
            "xlab": "log2 fold change",
            "ylab": "-log10 adjusted p-value",
            "title": f"DESeq2 results for {contrast[0]}",
            "data_labels": list(de_dict.keys()),
            "ymin": 0,
            "y_lines": [
                {"value": -np.log10(filter), "color": "#ff0000", "width": 1, "dash": "dash", "label": "significance"}
            ]
        },
        "data": reports
    }
    return mqc_report



def neighbors_report(adatas, spotlight=None, groups=None, ignore_self=True):
    self_type = "omitted" if ignore_self else "included"
    if spotlight:
        cell_types = spotlight    
    else:
        cell_types = {}
        for adata in adatas:
            if "cell_type" in adata.obs:
                ct_counts = adata.obs['cell_type'].value_counts()
                for ct, count in ct_counts.items():
                    if ct in cell_types:
                        cell_types[ct] += count
                    else:
                        cell_types[ct] = count
        cell_types = list(cell_types.keys())
    
    reports = []
    for cell_type in cell_types:
        sample_dict = {}
        for adata in adatas:
            # cell type interactions counts of immediate neighbors
            if "cell_type_interactions" in adata.uns:
                if cell_type not in adata.obs['cell_type'].cat.categories:
                    continue
                adata_cell_type_index = adata.obs['cell_type'].cat.categories.get_loc(cell_type)
                interactions = adata.uns["cell_type_interactions"][adata_cell_type_index]
                #construct dict of other cell types and interaction values
                other_cell_types = adata.obs['cell_type'].cat.categories.tolist()
                interaction_dict = {}
                for i, other_cell_type in enumerate(other_cell_types):
                    if ignore_self and other_cell_type == cell_type:
                        continue
                    interaction_dict[other_cell_type] = interactions[i]
                sample_dict[adata.obs['id'].unique()[0]] = interaction_dict
        sample_dict['cell_type'] = cell_type
        reports.append(sample_dict)

    # construct a df for export
    csv_report = pd.concat([pd.DataFrame(r).set_index('cell_type', append=True) for r in reports])
    csv_report.index.names = ['neighbor', 'cell_type']

    # recycle cell type as not needed in mqc report
    _cell_type = [r.pop('cell_type') for r in reports]

    mqc_report = {
        "id": "spatial_neighbors",
        "plot_type": "bar",
        "description": f"Cell type neighborhood across samples. The rate of self-neighborhood indicates clustering of a cell type. The rate of neighbors with other cell types (self omitted) indicates how these clusters interact with each other. Self is {self_type} in this report based on the pipeline settings.",
        "pconfig": {
            "title": "Cell type neighborhood across samples",
            "ylab": "Neighboring cell type share",
            "xlab": "Sample",
            "data_labels": cell_types
            },
        "data": reports
    }
    return mqc_report, csv_report

def diff_neighbors_report(neighbors_df:pd.DataFrame, group1=None, group2=None):
    # ALR transform using self neighborhood as the reference and pad for 0s
    report = neighbors_df.copy()
    for neighbor, cell_type in neighbors_df.index:
        self_value = neighbors_df.loc[(cell_type, cell_type)]
        alr = np.log((neighbors_df.loc[(neighbor, cell_type)] + 1e-10) / (self_value + 1e-10))
        report.loc[(neighbor, cell_type)] = alr
    
    # run ttest across samples for each neighbor-cell_type pair
    diff_neighbors = xsample_ttest(report, group1, group2)
    return diff_neighbors

def centrality_reports(adatas, spotlight=None, groups=None, uns_key='cell_type_centrality_scores'):
    #separately for each score as multiqc does not support data_labels for heatmaps
    memos = {
        'closeness_centrality': "This graph measure reflects how close the group is to other nodes.",
        'average_clustering': "This graph measure shows the degree to which nodes cluster together.",
        'degree_centrality': "This graph measure is the fraction of non-group members connected to group members"
    }

    adatas_scores = [a.uns[uns_key].keys().tolist() for a in adatas if uns_key in a.uns]
    if not adatas_scores:
        return {}
    scores = set.union(*[set(s) for s in adatas_scores])
    reports = {}
    for score in scores:
        sample_dict = {}
        for adata in adatas:
            if uns_key not in adata.uns or score not in adata.uns[uns_key]:
                continue
            centrality_scores = adata.uns[uns_key][score]
            available_cell_types = adata.obs['cell_type'].cat.categories.tolist()
            if spotlight:
                cell_types = [cell_type for cell_type in available_cell_types if cell_type in spotlight]
            else:
                cell_types = available_cell_types
            centrality_dict = {}
            for i, cell_type in enumerate(available_cell_types):
                if cell_type in cell_types:
                    centrality_dict[cell_type] = centrality_scores.iloc[i]
            sample_dict[adata.obs['id'].unique()[0]] = centrality_dict
            # differential scores if groups are provided
            csv_report = pd.DataFrame(sample_dict)
            if groups:
                try:
                    group1, group2 = groups
                    csv_report = xsample_ttest(csv_report, group1, group2)
                except Exception as e:
                    log.warning(f"Could not perform t-test for score {score}: {e}")

        mqc_report = {
            "id": f"{score}",
            "description": memos.get(score, f"Centrality score {score} across samples"),
            "plot_type": "table",
            "pconfig": {
                "ylab": "Sample",
                "xlab": "Cell type",
                "title": f"{score} centrality score across samples"
            },
            "data": sample_dict
        }
        reports[score] = [mqc_report, csv_report]
    return reports


def co_occurrence_report(adatas, spotlight=None, uns_key='cell_type_co_occurrence', summary='mean'):
    if spotlight:
        cell_types = spotlight    
    else:
        cell_types = {}
        for adata in adatas:
            if "cell_type" in adata.obs:
                ct_counts = adata.obs['cell_type'].value_counts()
                for ct, count in ct_counts.items():
                    if ct in cell_types:
                        cell_types[ct] += count
                    else:
                        cell_types[ct] = count
        cell_types = cell_types.keys()

    # check that number of cell types matches the co-occurrence matrix shape
    for adata in adatas:
        if len(adata.obs['cell_type'].cat.categories) != adata.uns[uns_key]['occ'].shape[0]:
            raise ValueError(f"Cell type categories do not match co-occurrence matrix dimensions for sample {adata.obs['id'].unique()[0]}")
    # check that the co-occurence intervals are same
    intervals = [adata.uns[uns_key]['interval'] for adata in adatas]
    if len(set(map(tuple, intervals))) > 1:
        log.warning(f"Co-occurrence intervals have mismatched dimensions: {intervals}")

    reports = []

    ct_dict = {}
    for ct in cell_types:
        # gather co-occurence for ct from each adata
        for adata in adatas:
            if uns_key not in adata.uns:
                continue
            if ct not in adata.obs['cell_type'].cat.categories:
                continue

            # gather all cell types in this adata
            adata_cell_types = adata.obs['cell_type'].cat.categories.tolist()
            # get the index of ct in current anndata
            ct_index = adata_cell_types.index(ct)
            # get its co-occurrence with other cell types
            co_occurrence = adata.uns[uns_key]['occ'][ct_index]
            # init a pandas dataframe with this cell and partner cell
            # types as composite index for further merging across samples
            sample = adata.obs['id'].unique()[0]
            # merge with previous samples on the cell type index, keeping all cell types seen in any sample
            co_df = pd.DataFrame(
                co_occurrence,
                index=pd.MultiIndex.from_tuples(
                    [(sample, ct, coct) for coct in adata_cell_types],
                    names=['sample', 'cell_type', 'co_cell_type'],
                ),
            )

            if ct in ct_dict:
                log.info(f"Merging co-occurrence data for cell type {ct} across samples.")
                ct_dict[ct] = pd.concat([ct_dict[ct], co_df])
            else:
                log.info(f"Initializing co-occurrence data for cell type {ct} with first sample.")
                ct_dict[ct] = co_df

    #return combined df if no summary requested, otherwise mqc
    csv_report = pd.concat([v for v in ct_dict.values()])
    
    for ct, df in ct_dict.items():
        df_summary = df.groupby(['co_cell_type']).agg(summary)
        report_dict = df_summary.to_dict(orient='index')
        # make list of points for each cell type, example {'stroma':[[0, 1.06],...[1, 0.95]]}
        report = {k:[[i,report_dict[k][i]] for i in report_dict[k].keys()] for k in report_dict.keys()}
        reports.append(report)

    mqc_report = {
        "id": "co_occurrence",
        "plot_type": "linegraph",
        "description": f"{summary} of Co-occurrence of cell types across samples. \
        Values greater than 1 indicate that the cell type is found \
        near other cell types more than average at a given distance.",
        "pconfig": {
            "title": "Cell type co-occurrence across samples",
            "ylab": "Occurring / expected",
            "xlab": "Interval",
            "data_labels": list(cell_types)
        },
        "data": reports
    }
    return mqc_report, csv_report



def plot_hist(df_pair, title=None, save=True):
    import matplotlib.pyplot as plt
    import seaborn as sns
    x = 8
    y = df_pair.shape[0] / 4

    plt.figure(figsize=(x, y))
    sns.heatmap(df_pair, cmap='coolwarm', fmt=".2f", linewidths=.5)
    plt.title(title)
    if save:
        plt.savefig(title.replace(" ", "_")+".png", dpi=300, bbox_inches='tight')

def xsample_ttest(df, group1, group2):
    res = df.copy()
    test = sp.stats.ttest_ind(res[group1], res[group2], axis=1)
    res['statistic'] = test.statistic
    res['pval'] = test.pvalue
    res.dropna(inplace=True)
    res['pval_adj'] = sp.stats.false_discovery_control(res['pval'], method='bh')
    res.sort_values('pval_adj', inplace=True)

    return res

def versions():
    with open ("versions.yml", "w") as f:
        f.write(f"{process}:\\n")
        f.write(f"    scipy: {sp.__version__}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
        f.write(f"    anndata: {ad.__version__}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    json: {json.__version__}\\n")


def get_vars(adatas, only=None):
    # find vars added by the staple pipeline
    vars = [x.uns['staple_meta_fields'].tolist() for x in adatas if 'staple_meta_fields' in x.uns]
    if (len(vars) == 0):
        log.warning("No added metadata fields found in any of the provided anndatas.")
        return None

    # drop id var from analysis vars
    for v in vars:
        if 'id' in v:
            v.remove('id')

    # find vars present in all samples
    common_vars = set(vars[0]).intersection(*vars[1:])
    log.info(f"Common added metadata fields across all samples: {common_vars}")
    if len(common_vars) == 0:
        log.warning("No common added metadata fields found across all provided anndatas.")

    res = common_vars

    if(only=='categorical'):
        # check that vars are categorical, same inside sample and differ across samples
        cats = {}
        for var in common_vars:
            is_categorical = all([isinstance(x.obs[var].dtype, pd.CategoricalDtype) for x in adatas])
            if not is_categorical:
                log.warning(f"Variable {var} is not categorical in all samples, skipping.")
                continue
            n_categories = [x.obs[var].nunique() for x in adatas]
            if len(set(n_categories)) != 1:
                log.warning(f"Variable {var} differs across individual samples .obs, skipping.")
                continue
            combined_categories = pd.concat([x.obs[[var,'id']] for x in adatas])
            levels = combined_categories.groupby(var, observed=True)['id'].unique()
            if (len(levels) != 2):
                log.warning(f"Can only contrast vars with 2 levels ({var} has {len(levels)}), skipping.")
                continue
            cats[var] = levels.to_dict()
        res = cats

    return res

def save_reports(mqc, res, name, mqc_reports="reports/mqc", reports="reports"):
    if mqc is not None:
        with open(f"{mqc_reports}/{name}_mqc.json","w") as f:
            json.dump(mqc, f, indent=4)
    if res is not None:
        res.to_csv(f"{reports}/{name}.csv")

if __name__ == '__main__':
    process       = "${task.process}"
    collected     = "${collected_items}"                 # these are whitespace separated paths to anndatas
    show          = int("${params.analyze.show_top}")    # how many top results to show
    cpus          = int("${task.cpus}")
    filter        = float("${params.analyze.filter}")    # p-value or adjusted p-value threshold for significance
    pb_vars       = "${params.analyze.pb_vars}"          # vars to use for pseudobulk grouping, comma-separated string
    only_spatial  = "${params.analyze.only_spatial}"\
                                    .lower() == 'true'   # only use spatially variable genes for ligrec and pseudobulk
    ignore_self   = "${params.analyze.ignore_self}"\
                                    .lower() == 'true'   # if true, ignore self interactions in neighbor analysis
    spotlight = "${params.analyze.spotlight}"            # a comma-separated string of cell type pairs to spotlight
    if spotlight:
        spotlight = [s.strip() for s in spotlight.split(',')]

    adata_paths = collected.split(" ")
    adatas = [ad.read_h5ad(path, backed="r") for path in adata_paths]

    #place all csv and mqc reports here
    reports_dir = "reports"
    os.makedirs(reports_dir, exist_ok=True)
    mqc_reports_dir = "reports/mqc"
    os.makedirs(mqc_reports_dir, exist_ok=True)

    # generate neighbors report
    try:
        log.info("Generating neighbors report.")
        neigh_mqc, neigh_csv = neighbors_report(adatas, spotlight=spotlight, ignore_self=ignore_self)
        save_reports(neigh_mqc, neigh_csv, "neighbors")
    except Exception as e:
        log.warning(f"Could not generate neighbors report: {e}")

    # generate centrality report - separately
    try:
        log.info("Generating centrality report.")
        centrality = centrality_reports(adatas, spotlight=spotlight)
        for score, report in centrality.items():
            save_reports(report[0], report[1], f"centrality_{score}")
    except Exception as e:
        log.warning(f"Could not generate centrality report: {e}")

    # print versions now because later may be too late
    versions()

    # generate reports using added meta vars
    all_vars = get_vars(adatas)
    log.info(f"All variables added upstream for analysis: {all_vars}")
    cats = get_vars(adatas, only='categorical')
    log.info(f"Cat variables suitable for cross-sample contrasts: {cats}")
    
    # prepare pseudobulk adata for deseq2, with all variables as grouping variables
    if pb_vars:
        # split comma-separated variables, trim whitespace, and preserve order without duplicates
        raw_vars = [v.strip() for v in pb_vars.split(",")]
        pb_vars = [v for v in dict.fromkeys(raw_vars) if v]
        log.info(f"Using specified variables for pseudobulk grouping: {pb_vars}")
    else:
        # start from all_vars as a list, preserving order and removing duplicates
        pb_vars = list(dict.fromkeys(all_vars)) if all_vars is not None else []
        log.info(f"No specific variables for pseudobulk grouping specified, using all variables: {pb_vars}")
        # ensure required grouping variables are present while preserving order
        if 'cell_type' not in pb_vars:
            pb_vars.append('cell_type')  # add cell type as grouping variable
        if 'id' not in pb_vars:
            pb_vars.append('id')
    pb_adata = pseudobulk_adatas(adatas, vars=pb_vars, only_spatial=only_spatial)
    if(only_spatial):
        log.info(f"Used only spatially variable genes for pseudobulk, resulting in {pb_adata.shape[1]} genes.")
        pb_adata.write(f"{reports_dir}/svg_pseudobulk.h5ad")
    else:
        log.info(f"Used all genes for pseudobulk, resulting in {pb_adata.shape[1]} genes.")
        pb_adata.write(f"{reports_dir}/pseudobulk.h5ad")
    
    # make reports
    supported_heatmap_tools = ['squidpy_ligrec', 'spacemarkers_LRscores', 'Moran_I']
    # if no cats found, just produce overall reports
    if cats is None or len(cats) == 0:
        # cycle through tools with supported heatmap reports
        for r in supported_heatmap_tools:
            try:
                res_mqc, res = heatmap_report(adatas, spotlight=spotlight, show=show, tool=r, filter=filter)
                save_reports(res_mqc, res, f"{r}_overall")
            except Exception as e:
                log.warning(f"Could not generate overall report for {r}: {e}")
        try:
            co_occ_mqc, co_occ_csv = co_occurrence_report(adatas, spotlight=spotlight)
            save_reports(co_occ_mqc, co_occ_csv, "co_occurrence_overall")
        except Exception as e:
            log.warning(f"Could not generate overall co-occurrence report: {e}")

    # for variables with 2 groups, perform contrast tests
    else:
        for var in cats.keys():
            groups = [x for x in cats[var]]
            group1 = cats[var][groups[0]].tolist()
            group2 = cats[var][groups[1]].tolist()
            
            #diff neighbors report
            try:
                log.info("Generating differential neighbors report.")
                res_mqc, res_csv = neighbors_report(adatas, spotlight=spotlight, ignore_self=False)
                diff_res = diff_neighbors_report(res_csv, group1=group1, group2=group2)
                save_reports(None, diff_res, f"neighbors_diff_{var}")
            except Exception as e:
                log.warning(f"Could not generate neighbors report for variable {var}: {e}")
            
            # cycle through tools with supported heatmap reports for each variable
            for r in supported_heatmap_tools:
                try:
                    res_mqc, res = heatmap_report(adatas, groups=[group1,group2],
                                                  spotlight=spotlight, show=show, tool=r, filter=filter)
                    save_reports(res_mqc, res, f"{r}_diff_{var}")
                except Exception as e:
                    log.warning(f"Could not generate {r} report for variable {var}: {e}")

            # differential centrality reports
            try:
                log.info("Generating differential centrality reports.")
                centrality_reports_res = centrality_reports(adatas, spotlight=spotlight, groups=[group1, group2])
                for score, report in centrality_reports_res.items():
                    save_reports(report[0], report[1], f"centrality_{score}_diff_{var}")
            except Exception as e:
                log.warning(f"Could not generate centrality report for variable {var}: {e}")

            # deseq2 on pseudobulks split by cell type
            try:
                contrasts = [var]+[k for k in cats[var].keys()]
                stratify_var = 'cell_type' if 'cell_type' in pb_adata.obs else None
                strata = pb_adata.obs[stratify_var].unique() if stratify_var else ['unstratified']
                de_results = {}
                for ct in strata:
                    mask = pb_adata.obs[stratify_var] == ct if stratify_var else np.array([True]*pb_adata.shape[0])
                    adata = pb_adata[mask].copy()
                    if adata.shape[0] < 2:
                        log.warning(f"Not enough samples for DESeq2 analysis for cell type {ct} with variable {var}, skipping.")
                        continue
                    de = pydeseq_results(adata, 
                                        spotlight=spotlight, 
                                        cpus=cpus, 
                                        design=f"~{var}", 
                                        contrast= contrasts)
                    de.summary()
                    deseq_res = de.results_df
                    deseq_res = deseq_res[deseq_res['padj'] <= filter]
                    if(deseq_res.empty):
                        log.warning(f"No significant DE genes found for {ct} cell type with variable {var}.")
                    else:
                        log.info(f"Deseq2 results for {ct} cell type with variable {var}:{deseq_res.shape[0]} significant genes found.")
                        de_results[ct] = deseq_res
                        deseq_res.to_csv(f"{reports_dir}/deseq2_diff_{var}_{ct}.csv")
                mqc_report = de_report(de_results, spotlight=spotlight, filter=filter, show=show, contrast=contrasts)
                save_reports(mqc_report, None, f"deseq2_diff_{var}",
                                mqc_reports_dir, reports_dir)
            except Exception as e:
                log.warning(f"Could not perform DESeq2 analysis for variable {var}: {e}")

            # co-occurence by groups (place each adata in the appropriate group based on its obs)
            try:
                co_occ_mqc, co_occ_csv = co_occurrence_report(adatas, spotlight=spotlight)
                co_occ_csv.insert(0, 'group', np.where(co_occ_csv.index.get_level_values('sample').isin(group1), groups[0], groups[1]))
                reports = []
                cell_types = co_occ_csv.index.get_level_values('cell_type').unique()
                for ct in cell_types:
                    ct_df = co_occ_csv.xs(ct, level='cell_type')
                    group1_df = ct_df[ct_df['group'] == groups[0]].drop(columns='group')
                    group2_df = ct_df[ct_df['group'] == groups[1]].drop(columns='group')
                    # compute median distance between groups for each co-occurring cell type
                    g1 = group1_df.groupby(['co_cell_type']).median()
                    g2 = group2_df.groupby(['co_cell_type']).median()
                    common_index = g1.index.intersection(g2.index)
                    g1, g2 = g1.loc[common_index], g2.loc[common_index]
                    z_diff = (g1 - g2)
                    report_dict = z_diff.to_dict(orient='index')
                    z_diff_report = {
                        k: [[i, report_dict[k][i]] for i in report_dict[k].keys()]
                        for k in report_dict.keys()
                    }
                    reports.append(z_diff_report)
                mqc_report = {
                    "id": f"co_occurrence_diff_{var}",
                    "plot_type": "linegraph",
                    "description": f"Median difference of co-occurrence across groups of variable {var}.",
                    "pconfig": {
                        "title": f"Co-occurrence difference by {var}",
                        "ylab": "Median difference",
                        "xlab": "Co-occurring cell type",
                        "data_labels": list(cell_types)
                    },
                    "data": reports
                }
                save_reports(mqc_report, co_occ_csv, f"co_occurrence_diff_{var}", mqc_reports_dir, reports_dir)
            except Exception as e:
                log.warning(f"Could not generate co-occurrence report for variable {var}: {e}")

    #wrapup
    for adata in adatas:
        adata.file.close()

