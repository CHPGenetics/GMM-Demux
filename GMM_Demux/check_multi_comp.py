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

