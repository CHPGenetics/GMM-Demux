from GMM_Demux import estimator
from GMM_Demux import classify_drops
import BitVector
import numpy as np

from sys import argv
import pandas as pd
from statistics import mean
from scipy.optimize import minimize
from scipy import optimize 
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint

# Returns the binary array representing all combinations of cells
def obtain_base_bv_array(sample_num):
    base_bv_array = []
    bv_array_segments = [0]

    for i in range(sample_num):
        base_bv_array.append(BitVector.BitVector(size = sample_num) )
        base_bv_array[-1][i] = 1

    bv_array_segments.append(len(base_bv_array))

    for i in range(1, sample_num):
        for j in range(bv_array_segments[-2], bv_array_segments[-1]):
            first_zero = sample_num - base_bv_array[j].reverse().next_set_bit()

            for k in range(first_zero, sample_num):
                base_bv_array.append(base_bv_array[j] | base_bv_array[k])
        bv_array_segments.append(len(base_bv_array))

    return [BitVector.BitVector(intVal = 0, size = sample_num)] + base_bv_array


def get_empty_bv(sample_num):
    return BitVector.BitVector(size = sample_num)


# Return true if the bit_pos th bit of bv is 1
def check_set_bit(bv, bit_pos, sample_num):
    mask = BitVector.BitVector(size = sample_num)
    mask[bit_pos] = 1
    
    return int(mask & bv) != 0


def init_mask(sample_num):
    mask = BitVector.BitVector(size = sample_num)
    return mask


def set_bit(bv, bit_pos):
    bv[bit_pos] = 1
    return bv


def gather_multiplet_rates(venn_values, SSM_rate_ary, sample_num):
    total_drops = 0
    total_singlets = 0
    total_MSMs = 0
    total_SSMs = 0

    for i in range(len(venn_values)):
        if i < sample_num:
            singlet_num = venn_values[i] * (1 - SSM_rate_ary[i])
            SSM_num = venn_values[i] - singlet_num

            total_SSMs += SSM_num
            total_singlets += singlet_num

        else:
            total_MSMs += venn_values[i]

        total_drops += venn_values[i]

    singlet_rate = total_singlets / total_drops
    SSM_rate = total_SSMs / total_drops
    MSM_rate = total_MSMs / total_drops

    return MSM_rate, SSM_rate, singlet_rate


def obtain_HTO_GEM_num(data_df, base_bv_array, sample_num):
    HTO_GEM_ary = []

    # Obtain hto numbers
    for i in range(1, len(base_bv_array)):
        #print(base_bv_array[i])
        GEM_num = 0

        for j in range(len(base_bv_array)):
            if (base_bv_array[i] & base_bv_array[j] == base_bv_array[i]):
                GEM_num += (data_df["Cluster_id"] == j).sum()

        HTO_GEM_ary.append(GEM_num)
    
    #print(HTO_GEM_ary)
    return HTO_GEM_ary


def experiment_params_wrapper(params, HTO_GEM_ary, sample_num, scaler, base_bv_array, operator):
    #print("***Params: ", params)
    scaled_params = params.copy()

    for i in range(len(scaler)):
        scaled_params[i] = operator(params[i], scaler[i])

    drop_num = round(scaled_params[0])
    capture_rate = scaled_params[1]
    cell_num_ary = [round(cell_num) for cell_num in scaled_params[2:]]

    return -estimator.compute_observation_probability(drop_num, capture_rate, cell_num_ary, HTO_GEM_ary, base_bv_array, sample_num)


def param_scaling(params, scaler, operator):
    for i in range(len(scaler)):
        params[i] = operator(params[i], scaler[i])

    return params


def compute_scaler(params):
    scaler = []

    # Bigger steps for the first 2 params
    scaler.append(800 / params[0])
    scaler.append(500 / params[1])

    sample_share = []
    for param in params[2:]:
        sample_share.append( param / sum(params[2:]) * 40000)

    for share, param in zip(sample_share, params[2:]):
        scaler.append(share / param)

    return scaler


def obtain_experiment_params(base_bv_array, HTO_GEM_ary, sample_num, estimated_total_cell_num, params0 = None):
    drop_num = 80000
    capture_rate = 0.5
    cell_num_ary = [estimated_total_cell_num / sample_num for i in range(sample_num)]

    if params0 is None:
        params0 = [drop_num, capture_rate, *cell_num_ary]

    lower_bound = [1, 0.0] + [0 for i in range(sample_num)]
    upper_bound = [np.inf, 1.0] + [np.inf for i in range(sample_num)]

    # Compute scaler
    scaler = compute_scaler(params0)

    #print("params before scaling:", params0)
    # Scale all values
    params0 = param_scaling(params0, scaler, lambda x, y: x * y)
    lower_bound = param_scaling(lower_bound, scaler, lambda x, y: x * y)
    upper_bound = param_scaling(upper_bound, scaler, lambda x, y: x * y)

    #print("params:", params0)
    #print("bounds:", lower_bound, upper_bound)

    bounds = Bounds(lower_bound, upper_bound)

    # Linear constraints are not scaled in current version.
    constraint_func = [
            {"type": "ineq", "fun": lambda x: sum([a/b for a,b in zip(x[2:], scaler[2:])]) - estimated_total_cell_num * 0.99},
            {"type": "ineq", "fun": lambda x: - sum([a/b for a,b in zip(x[2:], scaler[2:])]) + estimated_total_cell_num * 1.01}
            ]

    res = minimize(experiment_params_wrapper, params0, args=(HTO_GEM_ary, sample_num, scaler, base_bv_array, lambda x,y: x/y), method='SLSQP', constraints=constraint_func, bounds=bounds, options={'verbose': 1, 'eps': 10, 'ftol' : 0.0000001})

    #print(res)
    final_param = param_scaling(res.x, scaler, lambda x, y: x / y)

    #print(final_param)
    return final_param 


def obtain_HTO_cell_n_drop_num(data_df, base_bv_array, sample_num, estimated_total_cell_num, confidence_threshold):
    hto_num_ary = []
    HTO_num_ary = []
    estimated_drop_num_ary = []

    # Obtain hto numbers
    for i in range(sample_num):
        hto_num = 0

        for j in range(len(base_bv_array)):
            if check_set_bit(base_bv_array[j], i, sample_num):
                hto_num += (data_df["Cluster_id"] == j).sum()

        hto_num_ary.append(hto_num)

    for hto_a_idx in range(sample_num):
        for hto_b_idx in range(hto_a_idx + 1, sample_num):

            shared_num = 0

            for j in range(len(base_bv_array)):
                if check_set_bit(base_bv_array[j], hto_a_idx, sample_num) and check_set_bit(base_bv_array[j], hto_b_idx, sample_num):
                    shared_num += (data_df["Cluster_id"] == j).sum()

            if shared_num > 0:
                estimated_drop_num = estimator.drop_num_estimator(hto_num_ary[hto_a_idx], hto_num_ary[hto_b_idx], shared_num)
                estimated_drop_num_ary.append(estimated_drop_num)


    captured_drop_num = mean(estimated_drop_num_ary)
    #print("Drop num ary:", estimated_drop_num_ary)
    #print("Total drop num:", drop_num)

    diff_num = 10000000000

    #for capture_rate in np.arange(1.0, 0.0, -0.005):
        #print("test: ", capture_rate)

    for capture_rate in np.arange(1.0, 0.0, -0.005):

        HTO_num_ary = []

        for hto_idx in range(sample_num):
            hto_num = hto_num_ary[hto_idx]
            cell_num = estimator.cell_num_estimator(hto_num, captured_drop_num, capture_rate)
            HTO_num_ary.append(cell_num)

        total_cell_num = sum(HTO_num_ary)
        #print("total_cell_num: ", total_cell_num)
        #print("abs(total_cell_num - estimated_total_cell_num): ", abs(total_cell_num - estimated_total_cell_num))

        if (abs(total_cell_num - estimated_total_cell_num) > diff_num):
            break

        diff_num = abs(total_cell_num - estimated_total_cell_num)


    #print(HTO_num_ary)
    return HTO_num_ary, captured_drop_num / capture_rate, capture_rate

