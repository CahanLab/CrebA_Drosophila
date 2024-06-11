import sys 
import os 
from scipy import sparse, io
import pandas as pd 
import numpy as np 
import scanpy as sc
import pickle
import pySingleCellNet as pySCN

version_id = sys.argv[1]
os.chdir('../results/' + version_id + "/SCN_application")

data_types = [ name for name in os.listdir() if os.path.isdir(os.path.join('', name)) ]
[cgenesA, xpairs, tspRF] = pickle.load(open("../../../accessory_data/continuum_drosophila_embryonic_development/processed_data/SCN_classifier_obj.pickle", "rb"))

for data_type in data_types: 
    print(data_type)
    sparsematrix = io.mmread(data_type + '/raw_query_exp.txt')
    m_dense = sparsematrix.toarray()

    var_names = np.genfromtxt(data_type + '/raw_query_rownames.txt', dtype=str)
    col_names = np.genfromtxt(data_type + '/raw_query_colnames.txt', dtype=str)
    
    exp_df = pd.DataFrame(m_dense, columns=col_names, index=var_names)
    query_adata = sc.AnnData(exp_df.T)

    adVal = pySCN.scn_classify(query_adata, cgenesA, xpairs, tspRF, nrand = 0)
    class_matrix = adVal.to_df()
    SCN_label_mat = adVal.obs
    result = pd.concat([SCN_label_mat, class_matrix], axis=1)
    result.to_csv(data_type + "/SCN_classification.csv")


