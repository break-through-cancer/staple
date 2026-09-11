#!/usr/bin/env python
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


def make_one_adata(n=25, m=1000, pct_mito=0.1, sample_id='sample', with_metadata=True, na_pct=0.1, seed=42):
    np.random.seed(seed)

    # create anndata object n obs by m vars with random counts
    adata = ad.AnnData(X=np.random.poisson(1, (n, m)), 
                    obs=pd.DataFrame(index=[f'obs_{i}' for i in range(n)]), 
                    var=pd.DataFrame(index=[f'var_{j}' for j in range(m)]))
    # make sparse 
    adata.X = adata.X.astype(np.float32)

    # mark a percentage of the genes as mitochondrial by name prefix
    mito_genes = np.random.choice(adata.var_names, size=int(pct_mito*m), replace=False)
    adata.var_names = ['MT-' + name if name in mito_genes else name for name in adata.var_names]

    # add spatial dimensions to represent a 5x5 grid
    adata.obsm['spatial'] = np.array([[i // 5, i % 5] for i in range(n)])
    library_id = 'spatial_data'
    adata.uns['spatial'] = {
        library_id: {
            'scalefactors': {
                'tissue_hires_scalef': 1.0,
                'spot_diameter_fullres': 1.0
            },
            'images': {}
        }
    }
    sq.gr.spatial_neighbors(adata, spatial_key='spatial')


    # add an assignment with all outside cells being stroma surrounding
    # tumor core and a couple randomly assigned other of type other
    adata.obs['cell_type'] = 'cancer'
    outside_indices = [i for i in range(n) if adata.obsm['spatial'][i, 0] in [0, 4] or adata.obsm['spatial'][i, 1] in [0, 4]]
    for i in outside_indices:
        adata.obs.at[f'obs_{i}', 'cell_type'] = 'stroma'
    adata.obs.at[f'obs_{np.random.randint(0,n-1)}', 'cell_type'] = 'other'  # random cell
    adata.obs.at[f'obs_{np.random.randint(0,n-1)}', 'cell_type'] = 'other'  # another random cell
    
    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
    
    # attach cell type interaction report
    sq.gr.interaction_matrix(adata, cluster_key='cell_type')
    
    # add a random cell type assignment of 3 cell types: tumor, stroma, other
    cell_types = ['tumor', 'stroma', 'other']
    adata.obs['cell_type'] = np.random.choice(cell_types, size=n)
    adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'], categories=cell_types)

    # make some genes differentially expressed in stroma vs tumor for testing
    stroma = adata.obs['cell_type'] == 'stroma'
    tumor = adata.obs['cell_type'] == 'tumor'
    adata.X[stroma, :50] += np.random.poisson(5, (stroma.sum(), 50))
    adata.X[tumor, :50] += np.random.poisson(1, (tumor.sum(), 50))

    # set a small number of cell types to NA (as category since sq.gr needs that)
    na_sample_size = min(int(na_pct * adata.obs.shape[0]), adata.obs.shape[0])
    na_indices = np.random.choice(adata.obs.shape[0], size=na_sample_size, replace=False)
    adata.obs['cell_type'] = adata.obs['cell_type'].cat.add_categories('NA')
    adata.obs.loc[adata.obs.index[na_indices], 'cell_type'] = 'NA'

    # add sample id for testing
    adata.obs['id'] = sample_id

    # add a random response variable for testing (no variation within adata)
    adata.obs['response'] = np.random.choice(['responder', 'non-responder'])

    # add a continuous variable for testing
    adata.obs['age'] = np.random.randint(20, 100)

    # simulate staple behavior of added metadata from samplesheet
    if with_metadata:
        adata.uns['staple_meta_fields'] = ['response', 'id', 'age']

    # compute Moran's I (results stored in adata.uns['moranI']) for testing
    sq.gr.spatial_autocorr(adata, mode="moran", n_jobs=1)
    
    # mark some genes as spatially variable for testing
    adata.var['spatially_variable'] = adata.uns['moranI']['I'] > adata.uns['moranI']['I'].median()
    
    # rename spatial genes for easier debugging
    adata.var_names = ['spatial_' + g if s else g for g, s in zip(adata.var_names, adata.var['spatially_variable'])]
    
    # add dummy ligand-receptor interactions for testing. ligand receptor data
    # is adata.uns['ligrec_means'], adata.uns['ligrec_pvalues']
    # with gene-pairs in row index, and celltype pairs in columns
    # make some interactions from spatially variable genes, and some from 
    # non-spatially variable genes
    cell_types = adata.obs['cell_type'].cat.categories
    spatial = adata.var[adata.var['spatially_variable']].index
    non_spatial = adata.var[~adata.var['spatially_variable']].index
    gene_pairs = ['-'.join([g1, g2]) for g1, g2 in zip(spatial[:5], spatial[:5])] + \
        ['-'.join([g1, g2]) for g1, g2 in zip(non_spatial[:5], non_spatial[:5])]
    celltype_pairs = ['-'.join([ct1, ct2]) for ct1 in cell_types for ct2 in cell_types]
    ligrec_means = pd.DataFrame(np.random.rand(len(gene_pairs), len(celltype_pairs)), index=gene_pairs, columns=celltype_pairs)
    ligrec_pvalues = pd.DataFrame(np.random.rand(len(gene_pairs), len(celltype_pairs)), index=gene_pairs, columns=celltype_pairs)
    adata.uns['ligrec_means'] = ligrec_means
    adata.uns['ligrec_pvalues'] = ligrec_pvalues

    # recompute interactions using the final cell_type categories so that
    # the stored interaction matrix and centrality scores stay in sync
    sq.gr.interaction_matrix(adata, cluster_key='cell_type')

    # add co-occurence scores testing
    sq.gr.co_occurrence(adata, cluster_key='cell_type')

    # add centrality measures
    sq.gr.centrality_scores(adata, cluster_key='cell_type')

    return adata


def make_many_adata(num_adatas=2, n=25, m=1000, pct_mito=0.1, with_metadata=True, seeds=None):
    adatas = []
    if seeds is None:
        seeds = np.random.randint(0, 10000, size=num_adatas)  # different seed for each adata
    else:
        if len(seeds) != num_adatas:
            raise ValueError("Length of seeds must match num_adatas")
    for i in range(num_adatas):
        adata = make_one_adata(n=n, m=m, pct_mito=pct_mito, sample_id=f'sample_{i}',
                               with_metadata=with_metadata, na_pct=0.1, seed=seeds[i])
        adatas.append(adata)
    return adatas

if __name__ == "__main__":
    # nf params
    num_adatas = "${num_adatas}"
    with_metadata = "${with_metadata}".lower() == 'true'
    seeds = [s.strip() for s in "${seeds}".split(',')]
    if len(seeds) == 1 and seeds[0] == '':  # handle empty string case
        seeds = None
    else:
        seeds = [int(seed) for seed in seeds]
    #compose adatas and write to disk
    adatas = make_many_adata(num_adatas=int(num_adatas), n=25, m=1000,
                             pct_mito=0.1, with_metadata=with_metadata, seeds=seeds)
    for i, adata in enumerate(adatas):
        adata.write_h5ad(f'{i}_adata.h5ad')