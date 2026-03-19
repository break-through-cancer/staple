#!/usr/bin/env python
import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq


num_adatas = "${num_adatas}"

def make_one_adata(n=25, m=1000, pct_mito=0.1, sample_id='sample'):

    # create anndata object n obs by m vars with random counts
    adata = ad.AnnData(X=np.random.poisson(1, (n, m)), 
                    obs=pd.DataFrame(index=[f'obs_{i}' for i in range(n)]), 
                    var=pd.DataFrame(index=[f'var_{j}' for j in range(m)]))
    #make sparse 
    adata.X = adata.X.astype(np.float32)

    #mark a percentage of the genes as mitochondrial by name prefix
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

    # add sample id for testing
    adata.obs['id'] = sample_id
    # add a random response variable for testing
    adata.obs['response'] = np.random.choice(['responder', 'non-responder'])
    
    #simulate staple behavior of added metadata from samplesheet
    adata.uns['added_metadata_fields'] = ['response', 'id']
    
    # add Moran's I to obs for testing
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100, n_jobs=1)


    return adata


def make_many_adata(num_adatas=2, n=25, m=1000, pct_mito=0.1):
    adatas = []
    for i in range(num_adatas):
        adata = make_one_adata(n=n, m=m, pct_mito=pct_mito, sample_id=f'sample_{i}')
        adatas.append(adata)
    return adatas

if __name__ == "__main__":
    adatas = make_many_adata(num_adatas=int(num_adatas), n=25, m=1000, pct_mito=0.1)
    for i, adata in enumerate(adatas):
        adata.write_h5ad(f'{i}_adata.h5ad')