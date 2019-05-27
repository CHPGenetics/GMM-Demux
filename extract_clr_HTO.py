from scipy.io import mmread

from sys import argv
from sys import exit 

import pandas
import numpy as np
from scipy import stats

if len(argv) != 5:
    print("USAGE: script <dir_name> <num_HTO> <skip_num> <file_name>")
    exit()

dir_name = argv[1]

if dir_name[-1] != '/':
    dir_name = dir_name + '/'

cell_matrix = (mmread(dir_name + 'matrix.mtx'))
cell_matrix = cell_matrix.todense()

cell_names = open(dir_name + "barcodes.tsv").read().splitlines()

features = open(dir_name + "features.tsv").read().splitlines()

data = pandas.DataFrame(cell_matrix, features, cell_names).T

num_of_HTO = int(argv[2])
num_of_skip = int(argv[3])

HTO_data = data.loc[:,data.columns.values.tolist()[-(num_of_skip+num_of_HTO):-num_of_skip]]

HTO_columns = HTO_data.columns.to_series().astype(str)

HTO_columns_new = HTO_columns.apply(lambda x: x.split('\t')[0])

HTO_data.columns = HTO_columns_new

def normalize(array):
    array = array.astype(float).T
    for j in range(array.shape[0]):
        array_list = array[j].tolist()
        array_list = np.array([1.0 if float(i) < 0 else float(i) + 1.0 for i in array_list])
        geomean = stats.gmean(array_list)
        array[j] = np.log(np.true_divide(array_list, geomean))
    return array.T

norm_sample = normalize(HTO_data.values)

norm_df = pandas.DataFrame(norm_sample, index = HTO_data.index, columns = HTO_data.columns)

norm_df.to_csv(argv[4])

