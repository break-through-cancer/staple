#!/usr/bin/env python3
import json
import anndata as ad
import squidpy as sq
import logging
import scanpy as sc
import os
import pickle
from matplotlib import pyplot as plt
import numpy as np
from importlib.metadata import version

def default_ligrec(adata, par, **kwargs):
    if "gene_symbols" in kwargs:
        gene_symbols = kwargs.pop("gene_symbols")
    else:
        gene_symbols = None
    ligrec=sq.gr.ligrec(
        adata,
        n_perms=par["sq_gr_ligrec_nperms"],
        n_jobs=cpus,
        cluster_key="cell_type",
        copy=True,
        use_raw=False,
        corr_method="fdr_bh",
        interactions_params=par["sq_gr_ligrec_interactions_params"],
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        alpha=par["sq_gr_ligrec_alpha"],
        gene_symbols=gene_symbols,
        threshold=par["sq_gr_ligrec_threshold"])

    return ligrec

def default_ligrec_pl(ligrec, par, **kwargs):
    source_groups = kwargs.pop("source_groups")
    target_groups = kwargs.pop("target_groups")
    save = kwargs.pop("save")

    #reduce to sq_pl_ligrec_max_interactions for plotting
    n_top = par["sq_pl_ligrec_max_interactions"]-1
    sig_ligrec = ligrec["pvalues"][source_groups][target_groups].apply(min, axis=1) <= par["sq_pl_ligrec_pvalue"]
    sig_where = np.where(sig_ligrec)
    sig_means = ligrec["means"][source_groups][target_groups].iloc[sig_where].apply(max, axis=1)
    min_mean = sorted(sig_means, reverse=True)[n_top] if len(sig_means) > n_top else sig_means.min()
    ind = ligrec["means"][source_groups][target_groups].iloc[sig_where][sig_means>min_mean].index
    seed = par["seed"]

    repacked = ligrec.copy()
    for k in ["means", "pvalues", "metadata"]:
        repacked[k] = repacked[k].loc[ind]

    sq.pl.ligrec(repacked,
            source_groups=source_groups,
            target_groups=target_groups,
            save=save,
            title="", # cell names and title overlap sometimes
            alpha=0.01,
            swap_axes=False,
            seed=seed)
    plt.close('all')
    return None


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    log = logging.getLogger()

    adata_path = "${adata}"
    sample = "${prefix}"
    process = "${task.process}"
    cpus = int("${task.cpus}")


    par = {"sq_gr_ligrec_threshold": float("${params.sq_gr_ligrec_threshold}"),
        "sq_gr_ligrec_alpha": float("${params.sq_gr_ligrec_alpha}"),
        "sq_gr_ligrec_nperms": int("${params.sq_gr_ligrec_nperms}"),
        "sq_pl_ligrec_pvalue": float("${params.sq_pl_ligrec_pvalue}"),
        "sq_pl_ligrec_max_interactions": int("${params.sq_pl_ligrec_max_interactions}"),
        "seed": int("${params.seed}")
    }

    if '${params.sq_gr_ligrec_interactions_params}':
        par["sq_gr_ligrec_interactions_params"] = json.loads('${params.sq_gr_ligrec_interactions_params}')
    else:
        par["sq_gr_ligrec_interactions_params"] = None

    log.info(f"received params:{par}")
    log.info(f"loading {adata_path}")
    adata = ad.read_h5ad(adata_path)
    log.info(f"adata is {adata}")
    
    #get most abundant cell type from bayestme
    if 'cell_type' not in adata.obs.columns and 'bayestme_cell_type_counts' in adata.obsm:
        log.info("cell_type {} not found in adata.obs, calculating from bayestme_cell_type_counts")
        most_abundant = np.argmax(adata.obsm['bayestme_cell_type_counts'], axis=1)
        adata.obs['cell_type'] = most_abundant
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
        

    #squidpy insists on figure dir naming, not creating plot outdir as usually
    base_path = os.getcwd()
    save_path = os.path.join(sample, "ligrec")
    os.makedirs(save_path, exist_ok=True)
    os.chdir(save_path)

    #filter non-na cell_types
    adata = adata[~adata.obs["cell_type"].isna()].copy()
    clusters = adata.obs["cell_type"].unique()

    #normalize and log1p
    log.info("normalizing and log1p transforming data")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # try running ligrec
    res = None
    try:
        res = default_ligrec(adata, par)
    except Exception as e:
        log.error(f"default ligrec failed: {e}")
        try:
            res = default_ligrec(adata, par, gene_symbols="index")
        except Exception as e2:
            log.error(f"ligrec without gene_symbols failed: {e2}")

    if res is None:
        log.error("ligrec did not return results.")
    else:
        log.info("ligrec completed successfully, saving to pickle")
        # dictionary of pandas frames: means, pvalues, metadata
        pickle.dump(res, open("ligrec_interactions.pickle", "wb"))

        log.info("saving ligrec metadata")
        # just metadata, as it's easy to comprehend, one row per interaction
        res["metadata"].to_csv(
            "ligrec_metadata.csv"
        )

        log.info("saving ligrec interaction plot")
        for s in clusters:
            try:
                log.info(f"plotting ligrec for source cluster {s}")
                default_ligrec_pl(ligrec=res, source_groups=s, target_groups=clusters, par=par, save=f"source_{s}.png")
            except Exception as e:
                log.error(f"ligrec plot for source {s} failed: {e}")
                continue
        
        for t in clusters:
            try:
                log.info(f"plotting ligrec for target cluster {t}")
                default_ligrec_pl(ligrec=res, source_groups=clusters, target_groups=t, par=par, save=f"target_{t}.png")
            except Exception as e:
                log.error(f"ligrec plot for target {t} failed: {e}")
                continue

    os.chdir(base_path)
    with open ("versions.yml", "w") as f:
        f.write("{}:\\n".format(process))
        f.write("    squidpy: {}\\n".format(version("squidpy")))
        f.write("    anndata: {}\\n".format(version("anndata")))
        f.write("    scanpy: {}\\n".format(version("scanpy")))