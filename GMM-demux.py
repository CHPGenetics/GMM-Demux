# Just to suppress sklearn import imp warning
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import compute_venn
import param_estimator
import pandas as pd
import experiment_estimator
import classify_drops
from sys import argv
import sys

####### Parsing parameters and preparing data #######
data = pd.read_csv(argv[1], index_col = 0)
estimated_total_cell_num = int(argv[2])
confidence_threshold = float(argv[3])
store_path = "./result"

####### Run classifier #######
GEM_num = data.shape[0]
sample_num = data.shape[1]

base_bv_array = compute_venn.obtain_base_bv_array(sample_num)
print([int(i) for i in base_bv_array])
(high_array, low_array) = classify_drops.obtain_arrays(data)

# Obtain classification result
GMM_full_df, class_name_ary = \
        classify_drops.classify_drops(base_bv_array, high_array, low_array, sample_num, GEM_num, data.index, data.columns.values)

# Store results
classify_drops.store_full_classify_result(GMM_full_df, class_name_ary, store_path)
classify_drops.store_simplified_classify_result(GMM_full_df, class_name_ary, store_path, sample_num, confidence_threshold)

####### Estimate SSM #######
# Clean up bad drops
purified_df = classify_drops.purify_droplets(GMM_full_df, confidence_threshold)

# Infer parameters
(cell_num_ary, drop_num, capture_rate) = compute_venn.obtain_HTO_cell_n_drop_num(purified_df, base_bv_array, sample_num, estimated_total_cell_num, confidence_threshold)

SSM_rate_ary = [param_estimator.compute_SSM_rate_with_cell_num(cell_num_ary[i], drop_num) for i in range(sample_num)]

print("==============================Summary==============================")
rounded_cell_num_ary = [round(cell_num) for cell_num in cell_num_ary]
print("\n\nNumber of cells per subpopulation",rounded_cell_num_ary)
print("Total number of cells",sum(rounded_cell_num_ary))
print("Total number of drops",round(drop_num))
print("Capture Rate: %4.1f" % (capture_rate * 100) )


SSM_rate_df = pd.DataFrame(data=[SSM_rate_ary], columns = data.columns, index = ["Relative SSM rate per HTO"])
print("\n\n",SSM_rate_df)

MSM_rate, SSM_rate, singlet_rate = compute_venn.gather_multiplet_rates(list(bv_realcell_venn_dict.values()), SSM_rate_ary, sample_num)
summary_df = pd.DataFrame(data=[[MSM_rate, SSM_rate, singlet_rate]], columns = ["MSM","SSM","Singlet"], index = ["Overall Rates among all drops"])
print("\n\n",summary_df)

print("\n\nOverall relative SSM rate in non-MSM drops: %f" % experiment_estimator.compute_relative_SSM_rate(SSM_rate, singlet_rate) )

print("%i\t%5.2f\t%5.2f\t%5.2f\t%5.2f" % (drop_num, singlet_rate * 100, MSM_rate * 100, SSM_rate * 100, (SSM_rate / (SSM_rate + singlet_rate) ) * 100) , file=sys.stderr)
