# Just to suppress sklearn import imp warning
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import pandas as pd
from GMM_Demux import compute_venn
from GMM_Demux import estimator
from GMM_Demux import classify_drops
from GMM_Demux import GMM_IO
from sys import argv
import sys
import argparse
from tabulate import tabulate
from scipy.special import comb

def main():
    ####### Parsing parameters and preparing data #######
    parser = argparse.ArgumentParser(prog='GMM-demux', conflict_handler='resolve')

    # Positional arguments have * number of arguments atm.
    parser.add_argument('input_path', help = "The input path of mtx files from cellRanger pipeline.", nargs="*")
    parser.add_argument('hto_array', help = "Names of the HTO tags, separated by ','.", nargs="*")

    # Optional arguments.
    parser.add_argument("-k", "--skip", help="Load a full classification report and skip the mtx folder. Requires a path argument to the full report folder. When specified, the user no longer needs to provide the mtx folder.", type=str)
    parser.add_argument("-x", "--extract", help="Names of the HTO tag(s) to extract, separated by ','. Joint HTO samples are combined with '+', such as 'HTO_1+HTO_2'.", type=str)
    parser.add_argument("-o", "--output", help="The path for storing the Same-Sample-Droplets (SSDs). SSDs are stored in mtx format. Requires a path argument.", type=str, default="SSD_mtx")
    parser.add_argument("-f", "--full", help="Generate the full classification report. Requires a path argument.", type=str)
    parser.add_argument("-c", "--csv", help="Take input in csv format, instead of mmx format.", action='store_true')
    parser.add_argument("-t", "--threshold", help="Provide the confidence threshold value. Requires a float in (0,1). Default value: 0.8", type=float, default=0.8)
    parser.add_argument("-s", "--simplified", help="Generate the simplified classification report. Requires a path argument.", type=str)
    parser.add_argument("-u", "--summary", help = "Generate the statstic summary of the dataset. Including MSM, SSM rates. Requires an estimated total number of cells in the assay as input.", type=int)
    parser.add_argument("-r", "--report", help="Store the data summary report. Requires a file argument. Only executes if -u is set.", type=str)
    parser.add_argument("-e", "--examine", help="Provide the cell list. Requires a file argument. Only executes if -u is set.", type=str)
    parser.add_argument("-a", "--ambiguous", help="The estimated chance of having a phony GEM getting included in a pure type GEM cluster by the clustering algorithm. Requires a float in (0, 1). Default value: 0.05. Only executes if -e executes.", type=float, default=0.05)
    parser.add_argument("-rs", "--random_seed", help = "If provided, this is passed to the GaussianMixture random_state parameter.", default=None, type=int)
    
    print("==============================GMM-Demux Initialization==============================")

    args = parser.parse_args()

    confidence_threshold = args.threshold
    print("Confidence threshold:", confidence_threshold)


    # Classify droplets
    if not args.skip:
        # Overwrite the positional arguments
        parser.add_argument('input_path', help = "The input path of mtx files from cellRanger pipeline.")
        parser.add_argument('hto_array', help = "Names of the HTO tags, separated by ','.")
        args = parser.parse_args()

        input_path = args.input_path
        hto_array = args.hto_array.split(',')

        output_path = args.output
        print("Output directory:", output_path)

        #TODO: add CLR to csv data.
        if args.csv:
            full_df, GMM_df = GMM_IO.read_csv(input_path, hto_array)
        else:
            full_df, GMM_df = GMM_IO.read_cellranger(input_path, hto_array)
        
        GEM_num = GMM_df.shape[0]
        sample_num = GMM_df.shape[1]


        ####### Run classifier #######
        base_bv_array = compute_venn.obtain_base_bv_array(sample_num)
        #print([int(i) for i in base_bv_array])
        (high_array, low_array) = classify_drops.obtain_arrays(GMM_df, args.random_seed)

        # Obtain extract array.
        if args.extract:
            extract_id_ary = []
            tag_name_ary = []

            for tag_name in args.extract.split(','):
                tag_name_ary.append(tag_name.split('+') )

            for tag_ary in tag_name_ary:
                mask = compute_venn.init_mask(sample_num)
                for tag in tag_ary:
                    hto_idx = hto_array.index(tag)
                    bv = compute_venn.set_bit(mask, hto_idx)

                for idx in range(0, len(base_bv_array) ):
                    if base_bv_array[idx] == mask:
                        extract_id = idx

                extract_id_ary.append(extract_id)

        else:
            extract_id_ary = None 


        # Obtain classification result
        GMM_full_df, class_name_ary = \
                classify_drops.classify_drops(base_bv_array, high_array, low_array, sample_num, GEM_num, GMM_df.index, GMM_df.columns.values)

        # Store classification results
        if args.full:
            print("Full classification result is stored in", args.full)
            classify_drops.store_full_classify_result(GMM_full_df, class_name_ary, args.full)

        if args.simplified:
            ########## Paper Specific ############
            #purified_df = classify_drops.purify_droplets(GMM_full_df, confidence_threshold)
            ########## Paper Specific ############
            print("Simplified classification result is stored in", args.simplified)
            classify_drops.store_simplified_classify_result(GMM_full_df, class_name_ary, args.simplified, sample_num, confidence_threshold)
        
        # Clean up bad drops
        purified_df = classify_drops.purify_droplets(GMM_full_df, confidence_threshold)

        # Store SSD result
        print("MSM-free droplets are stored in folder", output_path, "\n")
        
        SSD_idx = classify_drops.obtain_SSD_list(purified_df, sample_num, extract_id_ary)
        SSD_df = GMM_IO.store_cellranger(full_df, SSD_idx, output_path)

        # Record sample names for summary report.
        sampe_names = GMM_df.columns

    # Parse the full report.
    else:
        GMM_full_df, sample_num, class_name_ary, sampe_names = classify_drops.read_full_classify_result(args.skip)
        base_bv_array = compute_venn.obtain_base_bv_array(sample_num)
        purified_df = classify_drops.purify_droplets(GMM_full_df, confidence_threshold)
        SSD_idx = classify_drops.obtain_SSD_list(purified_df, sample_num)


    ####### If extract is eanbled, other functions are disabled #######
    if args.extract:
        exit()


    ####### Estimate SSM #######
    if args.summary:
        # Count bad drops
        negative_num, unclear_num = classify_drops.count_bad_droplets(GMM_full_df, confidence_threshold)

        estimated_total_cell_num = args.summary

        # Infer parameters
        HTO_GEM_ary = compute_venn.obtain_HTO_GEM_num(purified_df, base_bv_array, sample_num)

        params0 = [80000, 0.5]

        for i in range(sample_num):
            params0.append(round(HTO_GEM_ary[i] * estimated_total_cell_num / sum(HTO_GEM_ary[:sample_num])))

        combination_counter = 0
        try:
            for i in range(1, sample_num + 1):
                combination_counter += comb(sample_num, i, True)
                HTO_GEM_ary_main = HTO_GEM_ary[0:combination_counter]
                params0 = compute_venn.obtain_experiment_params(base_bv_array, HTO_GEM_ary_main, sample_num, estimated_total_cell_num, params0)
        except:
            print("GMM cannot find a viable solution that satisfies the droplet formation model. SSM rate estimation terminated.")
            sys.exit(0)
                

        # Legacy parameter estimation
        #(cell_num_ary, drop_num, capture_rate) = compute_venn.obtain_HTO_cell_n_drop_num(purified_df, base_bv_array, sample_num, estimated_total_cell_num, confidence_threshold)
        (drop_num, capture_rate, *cell_num_ary) = params0

        SSM_rate_ary = [estimator.compute_SSM_rate_with_cell_num(cell_num_ary[i], drop_num) for i in range(sample_num)]
        rounded_cell_num_ary = [round(cell_num) for cell_num in cell_num_ary]
        SSD_count_ary = classify_drops.get_SSD_count_ary(purified_df, SSD_idx, sample_num)
        count_ary = classify_drops.count_by_class(purified_df, base_bv_array)
        MSM_rate, SSM_rate, singlet_rate = compute_venn.gather_multiplet_rates(count_ary, SSM_rate_ary, sample_num)

        # Generate report
        full_report_dict = {
            "#Drops": round(drop_num),
            "Capture rate": "%5.2f" % (capture_rate * 100),
            "#Cells": sum(rounded_cell_num_ary),
            "Singlet": "%5.2f" % (singlet_rate * 100),
            "MSM": "%5.2f" % (MSM_rate * 100),
            "SSM": "%5.2f" % (SSM_rate * 100),
            "RSSM": "%5.2f" % (estimator.compute_relative_SSM_rate(SSM_rate, singlet_rate) * 100),
            "Negative": "%5.2f" % (negative_num / GMM_full_df.shape[0] * 100),
            "Unclear": "%5.2f" % (unclear_num / GMM_full_df.shape[0] * 100)
            }
        full_report_columns = [
            "#Drops",
            "Capture rate",
            "#Cells",
            "Singlet",
            "MSM",
            "SSM",
            "RSSM",
            "Negative",
            "Unclear"
            ]

        full_report_df = pd.DataFrame(full_report_dict, index = ["Total"], columns=full_report_columns)

        print("==============================Full Report==============================")
        print(tabulate(full_report_df, headers='keys', tablefmt='psql'))
        print ("\n\n")
        print("==============================Per Sample Report==============================")
        sample_df = pd.DataFrame(data=[
            ["%d" % num for num in rounded_cell_num_ary],
            ["%d" % num for num in SSD_count_ary],
            ["%5.2f" % (num * 100) for num in SSM_rate_ary]
            ],
            columns = sampe_names, index = ["#Cells", "#SSDs", "RSSM"])
        print(tabulate(sample_df, headers='keys', tablefmt='psql'))

        if args.report:
            print("\n\n***Summary report is stored in folder", args.report)
            with open(args.report, "w") as report_file:
                report_file.write("==============================Full Report==============================\n")
            with open(args.report, "a") as report_file:
                report_file.write(tabulate(full_report_df, headers='keys', tablefmt='psql'))
            with open(args.report, "a") as report_file:
                report_file.write("\n\n")
                report_file.write("==============================Per Sample Report==============================\n")
            with open(args.report, "a") as report_file:
                report_file.write(tabulate(sample_df, headers='keys', tablefmt='psql'))


        # Verify cell type 
        if args.examine:
            print("\n\n==============================Verifying the GEM Cluster==============================")

            ambiguous_rate = args.ambiguous
            print("Ambiguous rate:", ambiguous_rate)

            simplified_df = classify_drops.store_simplified_classify_result(purified_df, class_name_ary, None, sample_num, confidence_threshold)

            cell_list_path = args.examine
            cell_list = [line.rstrip('\n') for line in open(args.examine)]
            cell_list = list(set(cell_list).intersection(simplified_df.index.tolist()))

            ########## Paper Specific ############
            #cell_list_df = pd.read_csv(args.examine, index_col = 0)
            #cell_list = cell_list_df.index.tolist()
            ########## Paper Specific ############

            MSM_list = classify_drops.obtain_MSM_list(simplified_df, sample_num, cell_list)

            GEM_num = len(cell_list)
            MSM_num = len(MSM_list)
            print("GEM count: ", GEM_num, " | MSM count: ", MSM_num)

            phony_test_pvalue = estimator.test_phony_hypothesis(MSM_num, GEM_num, rounded_cell_num_ary, capture_rate)
            pure_test_pvalue = estimator.test_pure_hypothesis(MSM_num, drop_num, GEM_num, rounded_cell_num_ary, capture_rate, ambiguous_rate)

            print("Phony-type testing. P-value: ", phony_test_pvalue)
            print("Pure-type testing. P-value: ", pure_test_pvalue)
            
            cluster_type = ""

            if phony_test_pvalue < 0.01 and pure_test_pvalue > 0.01:
                cluster_type = " pure"
            elif pure_test_pvalue < 0.01 and phony_test_pvalue > 0.01:
                cluster_type = " phony"
            else:
                cluster_type = "n unclear"

            print("Conclusion: The cluster is a" + cluster_type + " cluster.")

            ########## Paper Specific ############
            #estimated_phony_cluster_MSM_rate = estimator.phony_cluster_MSM_rate(rounded_cell_num_ary, cell_type_num = 2)
            #estimated_pure_cluster_MSM_rate = estimator.pure_cluster_MSM_rate(drop_num, GEM_num, rounded_cell_num_ary, capture_rate, ambiguous_rate)
            #print(str(estimated_phony_cluster_MSM_rate)+","+str(estimated_pure_cluster_MSM_rate)+","+str(MSM_num / float(GEM_num))+","+str(phony_test_pvalue)+","+str(pure_test_pvalue)+","+str(GEM_num / float(purified_df.shape[0]))+","+str(cluster_type), file=sys.stderr)
            ########## Paper Specific ############


if __name__ == "__main__":
    main()
