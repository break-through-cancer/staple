#!/usr/bin/env python3
# Cross-sample analysis template script
# Input is a space-delimited string with adatas

import os
import logging
import pandas as pd
import scipy as sp
import anndata as ad
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

    if pvalues is not None:
        # filter sample ligrec pairs by p-value
        sig_mask = (pvalues <= filter).any(axis=1)
        ligrecs = ligrecs[sig_mask]

    if groups is not None:
        ligrecs_ttest = xsample_ttest(ligrecs, groups[0], groups[1])
        ligrecs_ttest_sig = ligrecs_ttest[ligrecs_ttest['pval_adj'] <= filter]

        # if no t-test results found, fall back to mean across samples
        if(len(ligrecs_ttest_sig) == 0):
            log.warning(f"No significant differential interactions found for {tool}.")
            res = ligrecs_ttest.sort_values('pval', ascending=False)[samples]
            memo = f"Top {show} differential interactions across samples shown \
                as no significant (p_adj<={filter}) interactions were found \
                between {groups[0]} and {groups[1]}."
        else:
            res = ligrecs_ttest_sig.sort_values('pval')
            memo = f"Significant differential interactions (p_adj<={filter}) \
                between {groups[0]} and {groups[1]}. "
            pval_annot = ligrecs_ttest_sig['pval_adj'][:show]
            memo += f"There are {len(ligrecs_ttest_sig['pval_adj'])} significant \
                interactions (p_adj {min(pval_annot):.2e} to {max(pval_annot):.2e})."
            if len(ligrecs_ttest_sig['pval_adj']) > show:
                memo += f" Top {show} shown."
    else:
        ligrecs['mean'] = ligrecs.mean(axis=1)
        res = ligrecs.sort_values('mean', ascending=False)[samples]
        memo = f"Top {show} mean interactions across samples shown as no groups \
            were specified."
    
    #join multiindex of squidpy ligrec into single index
    if(isinstance(res.index, pd.MultiIndex)):
        res.index = ['-'.join(map(str, idx)) for idx in res.index]
    
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

def pseudobulk_from_adatas(adatas, filter=None):
    # collect pseudobulk expression for each cell type and each adata
    pseudobulk_dict = {}
    cell_types = set()
    cell_types = cell_types.union(*[set(adata.obs['cell_type'].cat.categories) for adata in adatas if 'cell_type' in adata.obs])
    for cell_type in cell_types:
        pseudobulk_dict[cell_type] = {}
        for adata in adatas:
            if 'cell_type' in adata.obs and cell_type in adata.obs['cell_type'].cat.categories:
                adata = adata.to_memory() if adata.isbacked else adata
                id = adata.obs['id'].unique()[0]
                subset = adata[adata.obs['cell_type'] == cell_type]
                bulk = subset.X.sum(axis=0).A1 if sp.sparse.issparse(subset.X) else subset.X.sum(axis=0)
                bulk = pd.Series(bulk, index=subset.var_names)
                pseudobulk_dict[cell_type][id] = bulk
                #check that raw integer counts are provided
                assert np.all(np.mod(bulk, 1) == 0), f"Non-integer counts found in {id} {cell_type} pseudobulk"
    return pseudobulk_dict

def pydeseq_results(pseudobulk_df, spotlight=None, groups=None, filter=0.05, cpus=1):
    inference = DefaultInference(n_cpus=cpus)
    # construct metadata dataframe for deseq2
    metadata = pd.DataFrame({
        'sample': pseudobulk_df.columns,
        'group': [1 if s in groups[0] else 2 for s in pseudobulk_df.columns]
    })
    metadata.set_index('sample', inplace=True)
    # run deseq2
    dds = DeseqDataSet(counts=pseudobulk_df.transpose(),
                       metadata=metadata,
                       design='~group',
                       inference=inference)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["group", 1, 2], inference=inference)
    #print summary
    ds.summary()
    res = ds.results_df

    # filter by p-value and sort
    if filter:
        res = res[res['padj'] <= filter]
    res.sort_values('padj', inplace=True)

    return res



def neighbors_report(adatas, spotlight=None):
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
    
    reports = []
    for cell_type in cell_types:
        sample_dict = {}
        for adata in adatas:
            if "cell_type_interactions" in adata.uns:
                if cell_type not in adata.obs['cell_type'].cat.categories:
                    continue
                adata_cell_type_index = adata.obs['cell_type'].cat.categories.get_loc(cell_type)
                interactions = adata.uns["cell_type_interactions"][adata_cell_type_index]
                #construct dict of other cell types and interaction values
                other_cell_types = adata.obs['cell_type'].cat.categories.tolist()
                interaction_dict = {}
                for i, other_cell_type in enumerate(other_cell_types):
                    interaction_dict[other_cell_type] = interactions[i]
                sample_dict[adata.obs['id'].unique()[0]] = interaction_dict
        reports.append(sample_dict)

    mqc_report = {
        "id": "spatial_neighbors",
        "plot_type": "bar",
        "description": "Cell type immediate neighborhood across samples. \
        The rate of self-neighborhood indicates clustering of a cell type. \
        The rate of neighbors with other cell types (self omitted) indicates \
            how these clusters interact with each other.",
        "pconfig": {
            "title": "Cell type neighborhood across samples",
            "ylab": "Neighboring cell type share",
            "xlab": "Sample",
            "data_labels": list(cell_types)
        },
        "data": reports
    }
    return mqc_report



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


def get_cat_vars(adatas):
    # find vars added by the staple pipeline
    vars = [x.uns['added_metadata_fields'].tolist() for x in adatas if 'added_metadata_fields' in x.uns]
    if (len(vars) == 0):
        log.warning("No added metadata fields found in any of the provided anndatas.")

    # drop id var from analysis vars
    for v in vars:
        if 'id' in v:
            v.remove('id')

    # find vars present in all samples
    common_vars = set(vars[0]).intersection(*vars[1:])
    log.info(f"Common added metadata fields across all samples: {common_vars}")
    if len(common_vars) == 0:
        log.warning("No common added metadata fields found across all provided anndatas.")

    
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

    return cats

if __name__ == '__main__':
    process = "${task.process}"
    collected = "${collected_items}"            # these are whitespace separated paths to anndatas
    show = int("${params.analyze.show_top}")    # how many top results to show
    cpus = int("${task.cpus}")
    filter = float("${params.analyze.filter}")  # p-value or adjusted p-value threshold for significance

    adata_paths = collected.split(" ")
    adatas = [ad.read_h5ad(path, backed="r") for path in adata_paths]
    spotlight = "${params.analyze.spotlight}" # this is a comma-separated string of cell type pairs to spotlight
    if ',' in spotlight:
        spotlight = spotlight.split(',')

    #place all mqc reports here
    reports_dir = "reports"
    os.makedirs(reports_dir, exist_ok=True)
    mqc_reports_dir = "reports/mqc"
    os.makedirs(mqc_reports_dir, exist_ok=True)

    # generate neighbors report
    try:
        neighbors = neighbors_report(adatas, spotlight=spotlight)
        with open(f"{mqc_reports_dir}/neighbors_mqc.json","w") as f:
            json.dump(neighbors, f, indent=4)
    except Exception as e:
        log.warning(f"Could not generate neighbors report: {e}")

    # print versions now because later may be too late
    versions()

    # generate reports using added meta vars
    cats = get_cat_vars(adatas)
    log.info(f"Variables and number of groups suitable for cross-sample analysis: {cats}")
    
    # make ligand-receptor reports
    # mqc report is for showing, but csv should have full data
    def save_reports(mqc, res, name, mqc_reports=mqc_reports_dir, reports=reports_dir):
        with open(f"{mqc_reports}/{name}_mqc.json","w") as f:
            json.dump(mqc, f, indent=4)
        res.to_csv(f"{reports}/{name}.csv")

    # if no cats found, just produce overall ligrec report
    if len(cats) == 0:
        try:
            res_mqc, res = heatmap_report(adatas, spotlight=spotlight, show=show, tool='squidpy_ligrec', filter=filter)
            save_reports(res_mqc, res, "ligrec_overall_mqc")
        except Exception as e:
            log.warning(f"Could not generate overall ligand-receptor report: {e}")
        try:
            lrs_mqc, lrs = heatmap_report(adatas, spotlight=spotlight, show=show, tool='spacemarkers_LRscores', filter=filter)
            save_reports(lrs_mqc, lrs, "lrscores_overall_mqc")
        except Exception as e:
            log.warning(f"Could not generate overall LR scores report: {e}")
        try:
            moran_mqc, moran = heatmap_report(adatas, spotlight=spotlight, show=show, tool='Moran_I', filter=filter)
            save_reports(moran_mqc, moran, "moranI_overall_mqc")
        except Exception as e:
            log.warning(f"Could not generate overall Moran's I report: {e}")
    else:
        # for variables with 2 groups, perform appropriate tests
        for var in cats.keys():
            groups = [x for x in cats[var]]
            group1 = cats[var][groups[0]].tolist()
            group2 = cats[var][groups[1]].tolist()
            
            #squidpy ligrec
            try:
                res_mqc, res = heatmap_report(adatas, groups=[group1,group2],
                                              spotlight=spotlight, show=show, tool="squidpy_ligrec", filter=filter)
                save_reports(res_mqc, res, f"ligrec_diff_{var}_results")
            except Exception as e:
                log.warning(f"Could not generate ligand-receptor report for variable {var}: {e}")
            # SpaceMarkers LR scores
            try:
                lrs_mqc, lrs = heatmap_report(adatas, groups=[group1,group2],
                                              spotlight=spotlight, show=show, tool='spacemarkers_LRscores', filter=filter)
                save_reports(lrs_mqc, lrs, f"lrscores_diff_{var}_results")
            except Exception as e:
                log.warning(f"Could not generate LR scores report for variable {var}: {e}")

            # Moran's I
            try:
                moran_mqc, moran = heatmap_report(adatas, groups=[group1,group2],
                                                  spotlight=spotlight, show=show, tool='Moran_I', filter=filter)
                save_reports(moran_mqc, moran, f"Moran_I_diff_{var}_results")
            except Exception as e:
                log.warning(f"Could not generate Moran's I report for variable {var}: {e}")

            # deseq2 on pseudobulks
            try:
                pseudobulks_by_ct = pseudobulk_from_adatas(adatas)
                for ct, pseudobulks in pseudobulks_by_ct.items():
                    pseudobulks = pd.DataFrame(pseudobulks)
                    pseudobulks = pseudobulks.fillna(0)  # fill missing values with 0
                    deseq_res = pydeseq_results(pseudobulks, groups=[group1,group2], filter=filter, cpus=cpus)
                    log.info(f"Deseq2 results for {ct} cell type with variable {var}: {deseq_res.shape[0]} significant genes found.")
                    if(deseq_res.empty):
                        log.warning(f"No significant DE genes found for {ct} cell type with variable {var}.")
                    else:
                        deseq_res.to_csv(f"{reports_dir}/deseq2_diff_{var}_results_{ct}.csv")
            except Exception as e:
                log.warning(f"Could not perform DESeq2 analysis for variable {var}: {e}")

    #wrapup
    for adata in adatas:
        adata.file.close()
