import classify_drops
import pandas as pd
from sys import argv
import numpy as np

def get_HTO_cell_idx(high_array, threshold):
    high_idx = np.argwhere(high_array > threshold)
    return high_idx


def compute_confidence(high_array, low_array, high_ary_idx, all_ary_idx):
    product = np.full(len(high_array[0]), 1.0)

    for idx in all_ary_idx:
        if idx in high_ary_idx:
            product = product * high_array[idx]
        else:
            product = product * low_array[idx]

    return product


def get_shared_cell_idx(high_array, low_array, high_ary_idx, all_ary_idx, threshold):
    product = compute_confidence(high_array, low_array, high_ary_idx, all_ary_idx)
    residual_idx = np.argwhere(product > threshold)
    return residual_idx


def get_HTO_cell_num(high_array, threshold):
    high_idx = get_HTO_cell_idx(high_array, threshold) 
    return len(high_idx)


def get_shared_cell_num(high_array, low_array, high_ary_idx, all_ary_idx, threshold):
    residual_idx = get_shared_cell_idx(high_array, low_array, high_ary_idx, all_ary_idx, threshold) 
    return len(residual_idx)


if __name__ == "__main__":

    data = pd.read_csv(argv[1], index_col = 0)

    (high_array, low_array) = classify_drops.obtain_arrays(data)

    multi_threshold = float(argv[2])
    single_threshold = float(argv[3])

    high_ary_idx = list(map(int, argv[4].split(',')))
    all_ary_idx = list(map(int, argv[5].split(',')))

    all_idx = []
    for idx in all_ary_idx:
        all_idx.append(get_HTO_cell_num(high_array[idx], single_threshold))

    print(all_idx)

    print(get_shared_cell_num(high_array, low_array, high_ary_idx, all_ary_idx, multi_threshold))
