from scipy.io import mmread, mmwrite
import pandas as pd
import os
from sys import argv
import gzip
from scipy import stats
import numpy as np
import os
from scipy.sparse import csr_matrix

pd.options.mode.chained_assignment = None
                                                                                                                        
def clr_norm(data_df):
    for hto in data_df.columns.values:
        compensated_values = data_df.loc[:,hto].values + 1
        gmean = stats.gmean(compensated_values)
        #print(gmean)

        data_df.loc[:,hto] = np.log(np.true_divide(compensated_values, gmean))

    return data_df


def read_csv(path, hto_array):
    full_df = pd.read_csv(path, index_col = 0)
    GMM_df = full_df.copy()
    GMM_df = GMM_df[hto_array]
    GMM_df = clr_norm(GMM_df)
    return full_df, GMM_df


def read_cellranger(path, hto_array):
    mtx_file = gzip.open(os.path.join(path, 'matrix.mtx.gz'), 'r')
    cell_matrix = (mmread(mtx_file))
    cell_matrix = cell_matrix.todense()

    name_file = gzip.open(os.path.join(path, 'barcodes.tsv.gz'), 'rt')
    cell_names = name_file.read().splitlines()
    cell_names = [name for name in cell_names]

    feature_file = gzip.open(os.path.join(path, 'features.tsv.gz'), 'rt')
    features = feature_file.read().splitlines()
    simple_features = [(feature.split('\t')[1] if (len(feature.split('\t')) > 1) else feature) for feature in features]

    full_df = pd.DataFrame(cell_matrix, simple_features, cell_names).T
    #print(full_df)

    data_df = full_df[hto_array]

    data_df = clr_norm(data_df)

    full_df.columns = features

    return full_df, data_df


def store_cellranger(data_df, SSD_idx, path):
    if not os.path.exists(path):
        os.makedirs(path)

    mtx_file = gzip.open(os.path.join(path, 'matrix.mtx.gz'), 'w')

    SSD_df = data_df.loc[SSD_idx,:]
    mmwrite(mtx_file, csr_matrix(SSD_df.T.values) )

    #print(SSD_df)

    feature_file = gzip.open(os.path.join(path, 'features.tsv.gz'), 'wt')
    for feature in SSD_df.columns.values:
        feature_file.write(feature + "\n")

    cell_file = gzip.open(os.path.join(path, 'barcodes.tsv.gz'), 'wt')
    for name in SSD_df.index.values:
        cell_file.write(name + "\n")
