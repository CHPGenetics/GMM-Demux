import compute_venn
import param_estimator
from sys import argv
from math import pow

def compute_multiplet_rates(cell_num, sample_num, drop_num):
    base_bv_array = compute_venn.obtain_base_bv_array(sample_num)
    HTO_num_ary = [round(cell_num / sample_num) for i in range(sample_num)]
    cell_mix_rate = [param_estimator.compute_mix_rate(drop_num, HTO_num) for HTO_num in HTO_num_ary]
    bv_simcell_venn_dict = compute_venn.compute_venn(HTO_num_ary, base_bv_array, cell_mix_rate, sample_num)
    (MSM_rate, SSM_rate, singlet_rate) = \
        compute_venn.gather_multiplet_rates(list(bv_simcell_venn_dict.values()), cell_mix_rate, sample_num)
    return MSM_rate, SSM_rate, singlet_rate


def compute_multiplet_rates_asymp(cell_num, sample_num, drop_num):
    no_drop_rate = (1 - 1 / drop_num)
    cells_per_sample = round(cell_num / sample_num)
    drop_with_cells = (1 - pow(no_drop_rate, cell_num)) * drop_num
    single_sample_drops = sample_num * (1 - pow(no_drop_rate, cells_per_sample) ) * drop_num * pow(no_drop_rate, (sample_num - 1) * cells_per_sample)
    singlet_drops = cell_num * pow(no_drop_rate, cell_num - 1)

    singlet_rate = singlet_drops / drop_with_cells
    single_sample_rate = single_sample_drops / drop_with_cells
    MSM_rate = 1 - single_sample_rate
    SSM_rate = single_sample_rate - singlet_rate
    return MSM_rate, SSM_rate, singlet_rate, drop_with_cells


def compute_relative_SSM_rate_asymp(cell_num, drop_num):
    no_drop_rate = (1 - 1 / drop_num)
    drop_with_cells = (1 - pow(no_drop_rate, cell_num)) * drop_num
    singlet_drops = cell_num * pow(no_drop_rate, cell_num - 1)

    return 1 - singlet_drops / drop_with_cells


def compute_relative_SSM_rate(SSM_rate, singlet_rate):
    return SSM_rate / singlet_rate


def get_min_hto_num(cell_num, drop_num, SSM_threshold, sample_num = 1):
    sample_num = 1
    while True:
        #(MSM_rate, SSM_rate, singlet_rate) = compute_multiplet_rates(cell_num, sample_num, drop_num)
        (MSM_rate, SSM_rate, singlet_rate, drop_with_cells) = compute_multiplet_rates_asymp(cell_num, sample_num, drop_num)
        relative_SSM_rate = compute_relative_SSM_rate(SSM_rate, singlet_rate)
        #print(relative_SSM_rate)
        if (relative_SSM_rate <= SSM_threshold):
            break;
        else:
            sample_num += 1

    return sample_num


if __name__ == "__main__":
    cell_num = float(argv[1])
    sample_num = int(argv[2])
    drop_num = float(argv[3])
    
    #print(compute_multiplet_rates(cell_num, sample_num, drop_num) )
    MSM_rate, SSM_rate, singlet_rate, drop_with_cells = compute_multiplet_rates_asymp(cell_num, sample_num, drop_num)

    no_drop_rate = (1 - 1 / drop_num)
    drop_with_cells = (1 - pow(no_drop_rate, cell_num)) * drop_num
    drop_per_sample = drop_with_cells * (1-MSM_rate) / sample_num

    print(SSM_rate / (SSM_rate + singlet_rate))
    print(param_estimator.compute_SSM_rate_with_cell_num(cell_num / sample_num, drop_num))

    print("singlet: ", singlet_rate)
    print("SSM: ", SSM_rate)
    print("MSM: ", MSM_rate)
    #print(compute_relative_SSM_rate(SSM_rate, singlet_rate))
    #print(compute_relative_SSM_rate_asymp(cell_num / sample_num, drop_num))
    
#    cell_num = float(argv[1])
#    drop_num = float(argv[2])
#    SSM_threshold = float(argv[3])

#    print(get_min_hto_num(cell_num, drop_num, SSM_threshold) )
