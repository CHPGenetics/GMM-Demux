import param_estimator
import check_multi_comp
import classify_drops
import BitVector
import numpy as np

from sys import argv
import pandas as pd
from statistics import mean
from matplotlib import pyplot as plt

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


# Return the Venn dictionary
def compute_venn(HTO_num_ary, base_bv_array, cell_mix_rate, sample_num):
    bv_simcell_dict = {}
    bv_simcell_venn_dict = {}

    for i in range(sample_num):
        bv_simcell_dict[int(base_bv_array[i])] = [HTO_num_ary[i] * (1 - cell_mix_rate[i])]

    #for bv in base_bv_array:
        #print(bv, " ", int(bv))

    for bv in base_bv_array[sample_num:]:
        bv_simcell_dict[int(bv)] = []
    #bv_simcell_dict = {int(bv) : [] for bv in base_bv_array}
    #print(bv_simcell_dict)
    #print(bv_simcell_dict)

    for bv in base_bv_array:
        bv_key = int(bv)
        for i in range(sample_num):
            i_bv = BitVector.BitVector(size = sample_num)
            i_bv[i] = 1
            intersect_bv = i_bv | bv # intersect here means Venn intersect not bv intersect
            if intersect_bv != bv:
                intersect_val = mean(bv_simcell_dict[bv_key]) * cell_mix_rate[i]
                intersect_bv_key = int(intersect_bv)
                bv_simcell_dict[intersect_bv_key].append(intersect_val)


    #print("Num of cells for each category:", bv_simcell_dict)


    for bv in base_bv_array:
        bv_key = int(bv)
        bv_simcell_venn_dict[bv_key] = mean(bv_simcell_dict[bv_key])

    # This is very important. It makes sure the inclusive relationship.
    for bv in base_bv_array:
        bv_key = int(bv)

        for i in range(sample_num):

            if not check_set_bit(bv, i, sample_num):
                bv_simcell_venn_dict[bv_key] *= 1 - cell_mix_rate[i]


    for bv in base_bv_array:
        bv_key = int(bv)
        bv_simcell_venn_dict[bv_key] = round(bv_simcell_venn_dict[bv_key])
       

    #print("Num of cells for Venn diagram:", bv_simcell_venn_dict)
    return bv_simcell_venn_dict


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
                estimated_drop_num = param_estimator.drop_num_estimator(hto_num_ary[hto_a_idx], hto_num_ary[hto_b_idx], shared_num)
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
            cell_num = param_estimator.cell_num_estimator(hto_num, captured_drop_num, capture_rate)
            HTO_num_ary.append(cell_num)

        total_cell_num = sum(HTO_num_ary)
        print("total_cell_num: ", total_cell_num)
        print("abs(total_cell_num - estimated_total_cell_num): ", abs(total_cell_num - estimated_total_cell_num))

        if (abs(total_cell_num - estimated_total_cell_num) > diff_num):
            break

        diff_num = abs(total_cell_num - estimated_total_cell_num)


    #print(HTO_num_ary)
    return HTO_num_ary, captured_drop_num / capture_rate, capture_rate

