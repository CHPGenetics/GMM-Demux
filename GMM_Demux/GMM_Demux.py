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

def main():
    ####### Parsing parameters and preparing data #######
    parser = argparse.ArgumentParser()
       
    parser.add_argument('input_path', help = "The input path of mtx files from cellRanger pipeline.")
    parser.add_argument('hto_array', help = "Names of the HTO tags, separated by ','.")
    parser.add_argument('cell_num', help = "Estimated total number of cells across all HTO samples.")

    parser.add_argument("-f", "--full", help="Generate the full classification report. Require a path argument.", type=str)
    parser.add_argument("-s", "--simplified", help="Generate the simplified classification report. Require a path argument.", type=str)
    parser.add_argument("-o", "--output", help="The path for storing the Same-Sample-Droplets (SSDs). SSDs are stored in mtx format. Require a path argument.", type=str)
    parser.add_argument("-r", "--report", help="Store the data summary report. Require a file argument.", type=str)

    args = parser.parse_args()

    input_path = args.input_path
    hto_array = args.hto_array.split(',')
    estimated_total_cell_num = int(args.cell_num)
    if args.output:
        output_path = args.output
    else:
        output_path = "GMM_Demux_mtx"


    full_df, GMM_df = GMM_IO.read_cellranger(input_path, hto_array)
    confidence_threshold = 0.8 

    ####### Run classifier #######
    GEM_num = GMM_df.shape[0]
    sample_num = GMM_df.shape[1]

    base_bv_array = compute_venn.obtain_base_bv_array(sample_num)
    #print([int(i) for i in base_bv_array])
    (high_array, low_array) = classify_drops.obtain_arrays(GMM_df)

    # Obtain classification result
    GMM_full_df, class_name_ary = \
            classify_drops.classify_drops(base_bv_array, high_array, low_array, sample_num, GEM_num, GMM_df.index, GMM_df.columns.values)

    # Store classification results
    if args.full:
        print("Full classification result is stored in", args.full)
        classify_drops.store_full_classify_result(GMM_full_df, class_name_ary, args.full)

    if args.simplified:
        print("Simplified classification result is stored in", args.simplified)
        classify_drops.store_simplified_classify_result(GMM_full_df, class_name_ary, args.simplified, sample_num, confidence_threshold)

    ####### Estimate SSM #######
    # Count bad drops
    negative_num, unclear_num = classify_drops.count_bad_droplets(GMM_full_df, confidence_threshold)
    # Clean up bad drops
    purified_df = classify_drops.purify_droplets(GMM_full_df, confidence_threshold)

    # Select target samples
    SSD_idx = classify_drops.obtain_SSD_list(purified_df, sample_num)

    # Store SSD result
    print("***MSM-free droplets are stored in folder", output_path, "\n")
    SSD_df = GMM_IO.store_cellranger(full_df, SSD_idx, output_path)

    # Infer parameters
    (cell_num_ary, drop_num, capture_rate) = compute_venn.obtain_HTO_cell_n_drop_num(purified_df, base_bv_array, sample_num, estimated_total_cell_num, confidence_threshold)

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
            "Unclear"]

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
        columns = GMM_df.columns, index = ["#Cells", "#SSDs", "RSSM"])
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


if __name__ == "__main__":
    main()
