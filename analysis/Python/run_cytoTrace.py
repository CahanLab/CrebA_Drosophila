import scanpy as sc
import numpy as np 
import pandas as pd
from scipy import sparse, io

import scvelo as scv
import cellrank as cr

from cellrank.tl.kernels import CytoTRACEKernel

import os 

import pickle

folders_list = os.listdir("../results/v18/wt13_cytoTrace")
folders_list = [x for x in folders_list if "DS_Store" not in x]

for temp_folder in folders_list:
    print(temp_folder)

    sparsematrix = io.mmread("../results/v18/wt13_cytoTrace/" + temp_folder + '/raw_query_exp.txt')
    m_dense = sparsematrix.toarray()
    var_names = np.genfromtxt("../results/v18/wt13_cytoTrace/" + temp_folder + '/raw_query_rownames.txt', dtype=str)
    col_names = np.genfromtxt("../results/v18/wt13_cytoTrace/" + temp_folder + '/raw_query_colnames.txt', dtype=str)
    exp_df = pd.DataFrame(m_dense, columns=col_names, index=var_names)
    query_adata = sc.AnnData(exp_df.T)

    query_meta_tab = pd.read_csv("../results/v18/wt13_cytoTrace/" + temp_folder + "/query_sample_tab.csv", index_col=0)
    query_adata.obs = query_meta_tab

    adata = query_adata.copy()
    sc.pp.filter_genes(adata, min_cells=10)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    # hvg annotation
    sc.pp.highly_variable_genes(adata)

    adata.layers["spliced"] = adata.X
    adata.layers["unspliced"] = adata.X

    n_neighbors = int(adata.shape[0] ** (1/2))
    scv.pp.moments(adata, n_neighbors=n_neighbors)
    ctk = CytoTRACEKernel(adata)

    #adata.write_h5ad("../results/v18/wt13_cytoTrace/" + temp_folder + "/ct_anndata.h5ad")
    pickle.dump(adata, open("../results/v18/wt13_cytoTrace/" + temp_folder + "/ct_anndata.pickle", "wb"))
    meta_tab = adata.obs
    meta_tab.to_csv("../results/v18/wt13_cytoTrace/" + temp_folder + "/cytotraced_sample_tab.csv")


