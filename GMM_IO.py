from scipy.io import mmread

from sys import argv
from sys import exit 

import pandas

if len(argv) != 7:
    print("USAGE: script <dir_name> <num_ADT> <num_HTO> <RNA_name> <ADT_name> <HTO_name>")
    exit()

dir_name = argv[1]

if dir_name[-1] != '/':
    dir_name = dir_name + '/'

cell_matrix = (mmread(dir_name + 'matrix.mtx'))
cell_matrix = cell_matrix.todense()

cell_names = open(dir_name + "barcodes.tsv").read().splitlines()

features = open(dir_name + "genes.tsv").read().splitlines()

data = pandas.DataFrame(cell_matrix, features, cell_names).T

num_of_ADT = int(argv[2])
num_of_HTO = int(argv[3])

#RNA_data = data.loc[:,data.columns.values.tolist()[:-(num_of_ADT+num_of_HTO)]]
#ADT_data = data.loc[:,data.columns.values.tolist()[-(num_of_ADT+num_of_HTO):-num_of_HTO]]
HTO_data = data.loc[:,data.columns.values.tolist()[-num_of_HTO:]]

#RNA_columns = RNA_data.columns.to_series().astype(str)
#ADT_columns = ADT_data.columns.to_series().astype(str)
HTO_columns = HTO_data.columns.to_series().astype(str)

#RNA_columns_new = RNA_columns.apply(lambda x: x.split('\t')[1])
#ADT_columns_new = ADT_columns.apply(lambda x: x.split('\t')[0])
HTO_columns_new = HTO_columns.apply(lambda x: x.split('\t')[0])

#RNA_data.columns = RNA_columns_new
#ADT_data.columns = ADT_columns_new
HTO_data.columns = HTO_columns_new

#RNA_data.to_csv(argv[4])
#ADT_data.to_csv(argv[5])
HTO_data.to_csv(argv[6])

#print(data.T)
