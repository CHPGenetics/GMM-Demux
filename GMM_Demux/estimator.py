from GMM_Demux import compute_venn
from math import pow
from math import log
from scipy.stats import binom
from scipy.stats import binom_test

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


def cell_num_estimator(a_num, captured_drop_num, capture_rate):
    estimated_drop_num = captured_drop_num / capture_rate
    base = 1 - 1 / estimated_drop_num
    power = 1 - a_num / captured_drop_num
    #print("power:", power)
    #print("base:", base)
    cell_num = log(power, base)
    return cell_num


def drop_num_estimator(a_num, b_num, shared_num):
    drop_num = a_num * b_num / shared_num 
    return drop_num 


def compute_shared_num(drop_num, A_num, B_num):
    A_rate = compute_mix_rate(drop_num, B_num) 
    #print("A_rate: ", A_rate)
    shared_num = A_rate *  A_num
    return shared_num


# Computes the rate of drops that have certain cells in them.
def compute_mix_rate(drop_num, cell_num):
    no_drop_rate = (1 - 1 / drop_num)
    cell_in_rate = 1 - pow(no_drop_rate, cell_num)
    return cell_in_rate


def compute_SSM_rate_with_cell_num(cell_num, drop_num):
    no_drop_rate = (1 - 1 / drop_num)
    singlet_drops = cell_num * pow(no_drop_rate, cell_num - 1)
    drop_with_cells = (1 - pow(no_drop_rate, cell_num)) * drop_num
    SSM_rate = 1 - singlet_drops / drop_with_cells
    return SSM_rate


def compute_SSD_num(drop_num, subject_cell_num, total_cell_num, ambiguous_rate = 0):
    no_drop_rate = (1 - 1 / drop_num)
    non_subject_num = total_cell_num - subject_cell_num
    SSD_prob = (1 - pow(no_drop_rate, subject_cell_num)) \
            * pow(no_drop_rate, int(non_subject_num * (1 - ambiguous_rate)))

    return int(SSD_prob * drop_num)


def compute_GEM_prob(drop_num, cell_num):
    return 1 - binom.pmf(0, cell_num, 1 / drop_num)


def phony_cluster_MSM_rate(cell_num_ary, cell_type_num = 2):
    assert(cell_type_num >= 2)

    total_cell_num = sum(cell_num_ary)
    sample_prob_ary = [cell_num / total_cell_num for cell_num in cell_num_ary]
    
    return 1 - sum([pow(sample_prob, cell_type_num) for sample_prob in sample_prob_ary])


def get_tau_cell_num(drop_num, total_cell_num, cluster_GEM_num, ambiguous_rate = 0.0):
    SSD_prob = cluster_GEM_num / drop_num 
    no_drop_rate = 1 - 1 / drop_num
    certain_rate = 1 - ambiguous_rate
    
    tau_cell_num = total_cell_num - int(log(SSD_prob \
            + pow(no_drop_rate, total_cell_num * certain_rate + cluster_GEM_num * ambiguous_rate) \
            , no_drop_rate) / certain_rate)

    return tau_cell_num


def pure_cluster_MSM_rate(drop_num, cluster_GEM_num, cell_num_ary, capture_rate, ambiguous_rate = 0):
    #print("==================")

    total_cell_num = sum(cell_num_ary)
    sample_prob_ary = [cell_num / total_cell_num for cell_num in cell_num_ary]

    cluster_GEM_num = cluster_GEM_num / capture_rate

    #print(cluster_GEM_num) 

    mix_cell_num = get_tau_cell_num(drop_num, total_cell_num, cluster_GEM_num)
    tau_cell_num = mix_cell_num 
    #tau_cell_num = (mix_cell_num - total_cell_num * ambiguous_rate) / (1 - ambiguous_rate)
    #tau_cell_num = get_tau_cell_num(drop_num, total_cell_num, cluster_GEM_num)

    #print("initial tau_cell_num: ", tau_cell_num)

    pure_type_GEM_num = compute_SSD_num(drop_num, tau_cell_num + (total_cell_num - tau_cell_num) * ambiguous_rate, total_cell_num, ambiguous_rate)
    while (pure_type_GEM_num > cluster_GEM_num):
        tau_cell_num -= 1
        pure_type_GEM_num = compute_SSD_num(drop_num, tau_cell_num + (total_cell_num - tau_cell_num) * ambiguous_rate, total_cell_num, ambiguous_rate)

    #print("final tau_cell_num: ", tau_cell_num)

    tau_num_ary = [int((tau_cell_num + (total_cell_num - tau_cell_num) * ambiguous_rate) * sample_prob) for sample_prob in sample_prob_ary]
    SSD_num_ary = [compute_SSD_num(drop_num, tau_num, total_cell_num) for tau_num in tau_num_ary]
    total_SSD_num = sum(SSD_num_ary)

    #print("tau GEM num: ", pure_type_GEM_num)
    #print("total num: ", compute_SSD_num(drop_num, total_cell_num, total_cell_num))

    return 1 - (total_SSD_num / pure_type_GEM_num)


def test_phony_hypothesis(cluster_MSM_num, cluster_GEM_num, cell_num_ary, capture_rate):
    MSM_rate = phony_cluster_MSM_rate(cell_num_ary)
    return binom_test(cluster_MSM_num / capture_rate, cluster_GEM_num / capture_rate, MSM_rate, "less")


def test_pure_hypothesis(cluster_MSM_num, drop_num, cluster_GEM_num, cell_num_ary, capture_rate, ambiguous_rate = 0):
    MSM_rate = pure_cluster_MSM_rate(drop_num, cluster_GEM_num, cell_num_ary, capture_rate, ambiguous_rate)
    print("Estimated MSM rate: ", MSM_rate)
    return binom_test(cluster_MSM_num / capture_rate, cluster_GEM_num / capture_rate, MSM_rate, "greater")


####Debuging Functions####
def debug_get_cell_num(drop_num, GEM_num, capture_rate):
    no_drop_rate = (1 - 1 / drop_num)
    return log(1 - GEM_num / (drop_num * capture_rate), no_drop_rate)


def debug_compute_doublet_num(drop_num, type_a_num, type_b_num):
    no_drop_rate = (1 - 1 / drop_num)
    SSD_prob = (1 - pow(no_drop_rate, type_a_num)) * (1 - pow(no_drop_rate, type_b_num))

    return int(SSD_prob * drop_num)


def debug_pure_cluster_MSM_rate(drop_num, tau_cell_num, sample_num_ary, capture_rate, ambiguous_rate = 0):
    total_cell_num = sum(sample_num_ary)
    sample_prob_ary = [sample_num / total_cell_num for sample_num in sample_num_ary]

    tau_num_ary = [int((tau_cell_num + (total_cell_num - tau_cell_num) * ambiguous_rate) * sample_prob) for sample_prob in sample_prob_ary]

    print(tau_num_ary)

    SSD_num_ary = [compute_SSD_num(drop_num, tau_num, total_cell_num) for tau_num in tau_num_ary]

    print(SSD_num_ary)

    total_SSD_num = sum(SSD_num_ary)

    print(total_SSD_num)

    pure_type_GEM_num = compute_SSD_num(drop_num, tau_cell_num + (total_cell_num - tau_cell_num) * ambiguous_rate, total_cell_num)
    #pure_type_GEM_num = compute_SSD_num(drop_num, tau_cell_num + (total_cell_num - tau_cell_num) * ambiguous_rate, total_cell_num * (1 - ambiguous_rate))

    print("Pure type num: ", int(pure_type_GEM_num * capture_rate))

    return 1 - (total_SSD_num / pure_type_GEM_num)
####End of Debuging Functions####


def compute_observation_probability(drop_num, capture_rate, cell_num_ary, HTO_GEM_ary, base_bv_array, sample_num):
    log_probability = 0

    GEM_prob_ary = []

    #print(drop_num, capture_rate, cell_num_ary, HTO_GEM_ary)
    
    for sample_idx in range(sample_num):
        ori_GEM_num = round(HTO_GEM_ary[sample_idx] / capture_rate)
        GEM_formation_prob = compute_GEM_prob(drop_num, cell_num_ary[sample_idx])
        GEM_prob_ary.append(GEM_formation_prob)

    for i in range(len(HTO_GEM_ary)):
        bv_idx = i + 1
        GEM_formation_prob = 1

        for sample_idx in range(sample_num):
            if compute_venn.check_set_bit(base_bv_array[bv_idx], sample_idx, sample_num):
                GEM_formation_prob *= GEM_prob_ary[sample_idx]

        ori_GEM_num = round(HTO_GEM_ary[i] / capture_rate)
        sample_binom_prob = binom.pmf(ori_GEM_num, drop_num, GEM_formation_prob)
        #print("***ori_GEM_num:", ori_GEM_num)
        #print("***GEM_formation_prob:", GEM_formation_prob)
        #print("***sample_binom_prob:", sample_binom_prob)
        #print("***log sample_binom_prob:", log(sample_binom_prob))
        #print("***log sample_binom_prob, corrected:", log(sample_binom_prob) * ((1/sample_num) ** (base_bv_array[bv_idx].count_bits() - 1)))
        log_probability += log(sample_binom_prob) * ((1/sample_num) ** (base_bv_array[bv_idx].count_bits() - 1))
        #probability *= sample_binom_prob

    return log_probability
