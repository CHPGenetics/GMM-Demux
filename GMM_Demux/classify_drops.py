import numpy as np

import pandas as pd
from scipy import stats
from sys import argv
from sklearn.mixture import GaussianMixture
import os
from math import log2

from GMM_Demux import check_multi_comp
from GMM_Demux import compute_venn


def obtain_arrays(data):
    gmm = []
    high_array = []
    low_array = []

    for i in range(data.shape[1]):
        X = data.iloc[:,i].values[:, np.newaxis]

        # GMM values
        gmm.append(GaussianMixture(2).fit(X))
        x = np.linspace(-6, 6, 1000)[:, np.newaxis]
        logprob= gmm[-1].score_samples(x)
        responsibilities = gmm[-1].predict_proba(x)
        pdf = np.exp(logprob)
        pdf_individual = responsibilities * pdf[:, np.newaxis]
        #print(pdf_individual)

        #print(gmm[-1].means_)

        # Extract prob
        high_idx = np.argmax(gmm[-1].means_, axis=0)[0]
        post_prob = gmm[-1].predict_proba(X)
        high_array.append(post_prob[np.arange(post_prob.shape[0]), np.full(post_prob.shape[0], high_idx)])
        low_array.append(np.full(post_prob.shape[0], 1.0) - high_array[-1])

    return high_array, low_array


def classify_drops(base_bv_array, high_array, low_array, sample_num, GEM_num, index, column):

    classified_ary = np.full(GEM_num, 0)
    confidence_ary = np.full(GEM_num, 0.0)

    class_name_array = ["negative"]
    all_idx_ary = []

    for i in range(sample_num):
        all_idx_ary.append(i)

    # Detailed classification
    for i in range(len(base_bv_array)):
        bv = base_bv_array[i]
        high_idx_ary = []
        name = ""

        for j in range(sample_num):
            if compute_venn.check_set_bit(bv, j, sample_num):
                high_idx_ary.append(j)
                name += (column[j] + "-")
        
        if name != "":
            name = name[:-1]
            class_name_array.append(name)

        tmp_confidence_ary = check_multi_comp.compute_confidence(high_array, low_array, high_idx_ary, all_idx_ary)
        update_idx = (tmp_confidence_ary > confidence_ary).nonzero()[0]
        confidence_ary[update_idx] = tmp_confidence_ary[update_idx]
        classified_ary[update_idx] = i


    GMM_array = np.column_stack(tuple([classified_ary, confidence_ary]))
    GMM_full_df = pd.DataFrame(data=GMM_array, index = index, columns = ["Cluster_id", "Confidence"])
    GMM_full_df["Cluster_id"] = GMM_full_df["Cluster_id"].astype(int)

    return GMM_full_df, class_name_array


def read_full_classify_result(path):
    print(path)
    classify_file_name = os.path.join(path, "GMM_full.csv")
    config_file_name = os.path.join(path, "GMM_full.config")

    full_df = pd.read_csv(classify_file_name, index_col = 0)
    config_df = pd.read_csv(config_file_name, header=None)

    sample_num = int(log2(config_df.shape[0]))

    return full_df, sample_num, config_df.index.tolist(), [config_df.index[i] for i in range(sample_num)]


# Store full classification result
def store_full_classify_result(data_df, class_name_array, path):

    if not os.path.exists(path):
        os.makedirs(path)

    classify_file_name = os.path.join(path, "GMM_full.csv")
    config_file_name = os.path.join(path, "GMM_full.config")

    data_df.to_csv(classify_file_name)

    with open(config_file_name, 'w') as f:
        for i in range(len(class_name_array)):
            f.write("%s, %s\n" % (i, class_name_array[i]))


# Store simplified classification result
def store_simplified_classify_result(data_df, class_name_array, path, sample_num, confidence_threshold):

    simplified_df = data_df.copy()
    #print(simplified_df)
    MSM_idx = data_df.index[(data_df["Cluster_id"] > sample_num).nonzero()[0]]
    #print(MSM_idx)
    simplified_df.loc[MSM_idx, "Cluster_id"] = sample_num + 1
    unclear_idx = data_df.index[(data_df["Confidence"] < confidence_threshold).nonzero()[0]]
    #print(unclear_idx)
    simplified_df.loc[unclear_idx, "Cluster_id"] = sample_num + 2
    #print(simplified_df)

    if path != None:
        if not os.path.exists(path):
            os.makedirs(path)

        classify_file_name = os.path.join(path, "GMM_simplified.csv")
        config_file_name = os.path.join(path, "GMM_simplified.config")

        simplified_df.to_csv(classify_file_name)

        simplified_name_array = [class_name_array[i] for i in range(sample_num + 1)]
        simplified_name_array.append("MSM")
        simplified_name_array.append("Unclear")

        with open(config_file_name, 'w') as f:
            for i in range(len(simplified_name_array)):
                f.write("%s, %s\n" % (i, simplified_name_array[i]))

    return simplified_df


def purify_droplets(data_df, confidence_threshold):
    drop_idx = data_df.index[((data_df["Confidence"] < confidence_threshold) | (data_df["Cluster_id"] == 0)).to_numpy().nonzero()[0]]
    purified_df = data_df.drop(drop_idx)
    return purified_df


def count_bad_droplets(data_df, confidence_threshold):
    negative_num = (data_df["Cluster_id"] == 0).sum()
    unclear_num = (data_df["Confidence"] < confidence_threshold).sum()
    return negative_num, unclear_num


def obtain_SSD_list(data_df, sample_num, class_id_ary = None):
    if class_id_ary is not None:
        SSD_idx = []
        for class_id in class_id_ary:
            SSD_idx.extend(data_df.index[data_df["Cluster_id"] == class_id])
    else:
        SSD_idx = data_df.index[data_df["Cluster_id"] <= sample_num]

    return SSD_idx


def obtain_MSM_list(data_df, sample_num, idx_list = None):
    if idx_list == None:
        MSM_idx = data_df.index[data_df["Cluster_id"] == sample_num + 1]
    else:
        selected_df = data_df.loc[idx_list]
        MSM_idx = selected_df.index[selected_df["Cluster_id"] == sample_num + 1]

    return MSM_idx


def count_by_class(data_df, base_bv_array):
    count_ary = []
    for i in range(1, len(base_bv_array) + 1):
        count_ary.append((data_df["Cluster_id"] == i).sum())

    return count_ary


def get_SSD_count_ary(data_df, SSD_idx, sample_num):
    SSD_df = data_df.loc[SSD_idx,:]
    SSD_count_ary = []

    for i in range(1, sample_num + 1):
        SSD_count_ary.append((SSD_df["Cluster_id"] == i).sum())

    return SSD_count_ary


