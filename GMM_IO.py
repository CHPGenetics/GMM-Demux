from scipy.io import mmread
import pandas as pd
import os
from sys import argv
import tarfile
import gzip
from scipy import stats
import numpy as np

pd.options.mode.chained_assignment = None
                                                                                                                        
def clr_norm(data_df):
    for hto in data_df.columns.values:
        compensated_values = data_df.loc[:,hto].values + 1
        gmean = stats.gmean(compensated_values)
        print(gmean)

        data_df.loc[:,hto] = np.log(np.true_divide(compensated_values, gmean))

    return data_df


def read_cellranger(path, hto_array):
    mtx_file = gzip.open(os.path.join(path, 'matrix.mtx.gz'), 'r')
    cell_matrix = (mmread(mtx_file))
    cell_matrix = cell_matrix.todense()

    name_file = gzip.open(os.path.join(path, 'barcodes.tsv.gz'), 'r')
    cell_names = name_file.read().splitlines()
    cell_names = [name.decode("utf-8") for name in cell_names]

    feature_file = gzip.open(os.path.join(path, 'features.tsv.gz'), 'r')
    features = feature_file.read().splitlines()
    features = [feature.decode("utf-8").split('\t')[1] for feature in features]

    full_df = pd.DataFrame(cell_matrix, features, cell_names).T
    #print(full_df)

    data_df = full_df[hto_array]

    data_df.to_csv("debug.csv")

    data_df = clr_norm(data_df)

    print(data_df)

    return full_df, data_df


dir_name = argv[1]
hto_array = argv[2].split(',')
print(hto_array)

read_cellranger(dir_name, hto_array)
