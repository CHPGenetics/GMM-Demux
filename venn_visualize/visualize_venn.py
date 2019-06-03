from GMM_Demux import estimator
from GMM_Demux.compute_venn import obtain_base_bv_array
from GMM_Demux.compute_venn import check_set_bit 
import pyvenn.venn
import BitVector
import numpy as np

from sys import argv
import pandas as pd
from statistics import mean
from matplotlib import pyplot as plt


# Compute the Venn dictionary
def compute_venn(HTO_num_ary, base_bv_array, cell_mix_rate, sample_num):
    bv_simcell_dict = {}
    bv_simcell_venn_dict = {}

    for i in range(sample_num):
        bv_simcell_dict[int(base_bv_array[i])] = [HTO_num_ary[i] * (1 - cell_mix_rate[i])]
        print(int(base_bv_array[i]))

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
    print(bv_simcell_venn_dict)
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


if __name__ == "__main__":

    def run_venn(venn_num, figsize, title, *args):

        method_to_call = getattr(pyvenn.venn, "venn" + str(venn_num))
        #print(args)

        plt.ioff()
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        fig, ax = method_to_call(*args, figsize=figsize, fontsize=14)
        plt.title(title, fontsize=17)
        plt.show()


    # Fill in the dummy array for Venn Diagram
    def plot_venn(bv_ary, bv_simcell_dict, venn_num, title=""):
        venn_dummy_array = []
        venn_name_array = []
        for i in range(venn_num):
            venn_dummy_array.append([])
            venn_name_array.append("hto" + str(i+1))

        dummy_num = 0
        for bv in bv_ary:
            bv_key = int(bv)
            for i in range(bv_simcell_dict[bv_key]):

                for i in range(venn_num):

                    if check_set_bit(bv, i, sample_num):
                        venn_dummy_array[i].append(dummy_num)
               
                dummy_num += 1


        labels = pyvenn.venn.get_labels(venn_dummy_array, ['number', 'logic'])
        for key in labels.keys():
            labels[key] = labels[key].split(": ")[1]
        print(labels)
        run_venn(sample_num, (9,7), title, labels, venn_name_array)


    def dump_to_csv(base_bv_array, bv_realcell_venn_dict, bv_simcell_venn_dict, sample_num):
        entries = []

        for key in base_bv_array:
            name = ""
            for i in range(sample_num):
                if check_set_bit(key, i, sample_num):
                    if name != "":
                        name += "-"
                    
                    name += ("hto" + str(i + 1))

            entries.append((name, bv_realcell_venn_dict[int(key)], bv_simcell_venn_dict[int(key)]))

        labels = ["sample combination", "observed counts", "model-derived counts"]

        df = pd.DataFrame.from_records(entries, columns = labels)

        df.to_csv("venn.csv", index=False)
            

    GMM_df = pd.read_csv(argv[1], index_col = 0)
    sample_num = int(argv[2])
    params = argv[3].split(',')
    confidence_threshold = float(argv[4])

    assert(len(params) == sample_num + 2)

    sim_drop_num = round(float(params[0]))
    sim_capture_rate = float(params[1])
    GEM_cell_ary = []
    for param in params[2:]:
        GEM_cell_ary.append(round(float(param) ) )

    base_bv_array = obtain_base_bv_array(sample_num)[1:]

    GEM_prob_ary = [estimator.compute_GEM_prob(sim_drop_num, cell_num) for cell_num in GEM_cell_ary]

    # Compute the Venn dictionary
    bv_simcell_venn_dict = compute_venn(GEM_cell_ary, base_bv_array, GEM_prob_ary, sample_num)
    print(bv_simcell_venn_dict)

    # Account capture rate
    for bv in bv_simcell_venn_dict:
        bv_simcell_venn_dict[bv] = int(round(bv_simcell_venn_dict[bv] * sim_capture_rate))



    # Generate the bv_simcell_dict from a real sample
    bv_realcell_venn_dict = {}
    all_idx_ary = [i for i in range(sample_num)]
    for i in range(0, len(base_bv_array)):
        bv_key = int(base_bv_array[i])
        #bv_realcell_venn_dict[bv_key] = ((GMM_df["Cluster_id"] == i + 1) & (GMM_df["Confidence"] > confidence_threshold)).sum()
        bv_realcell_venn_dict[bv_key] = (GMM_df["Cluster_id"] == i + 1).sum()
    

    plot_venn(base_bv_array, bv_realcell_venn_dict, sample_num, "GMM-Demux-Observed Venn Diagram")
    print(bv_simcell_venn_dict)
    plot_venn(base_bv_array, bv_simcell_venn_dict, sample_num, "Model-Derived Venn Diagram")
    #dump_to_csv(base_bv_array, bv_realcell_venn_dict, bv_simcell_venn_dict, sample_num)
    #SSM_rate_ary = [param_estimator.compute_SSM_rate_with_cell_num(HTO_num_ary[i], drop_num) for i in range(sample_num)]
    #print(gather_multiplet_rates(list(bv_simcell_venn_dict.values()), SSM_rate_ary, sample_num))
