#!/usr/bin/env python3

import argparse
import sys
import re
import os
from glob import glob
from collections import OrderedDict
import json


def check_arg(args=None):
    parser = argparse.ArgumentParser(
        prog="kmerfinder_hits.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="extract the best and second hits found in kmerfinder ",
    )
    parser.add_argument(
        "-p",
        "--path",
        required=True,
        help="Insert path of results.txt file like /home/user/Service_folder/ANALYSIS/07-kmerfinder",
    )

    parser.add_argument(
        "-o", "--output", required=True, help="The output folder to save  hits"
    )

    # Example: python3 parse_kmerfinder.py -p /path_to_kmer_out_results -o p_kmer_hits
    return parser.parse_args()


def hits_look_up(file_txt):
    """Function to extract the best and the second hits from kmer result"""
    hits = {}
    with open(file_txt) as fh:
        lines = fh.readlines()
    parameters = lines[0].strip().split("\t")
    # check if get at least one best hits
    try:
        best_hit_data = lines[1].strip().split("\t")
    except IndexError:
        return hits
    hits["best_hit"] = OrderedDict()
    # get the second hits if exists
    try:
        second_hit_data = lines[2].strip().split("\t")
        hits["second_hit"] = OrderedDict()
    except IndexError:
        pass

    for idx in range(len(parameters)):
        hits["best_hit"][parameters[idx]] = best_hit_data[idx]
        if "second_hit" in hits:
            hits["second_hit"][parameters[idx]] = second_hit_data[idx]

    return hits


def save_data_to_json(data, file_name):
    with open(file_name, "w") as fh:
        json.dump(data, fh)
    return


if __name__ == "__main__":

    version = "v 0.1.0."
    args = check_arg(sys.argv[1:])

    path = args.path
    out_folder = args.output
    sample_list = glob(os.path.join(path, "*_results.txt"))

    kmer_all = {}
    for s_file in sample_list:
        tmp_s = re.search(r"(\w+)_results\.txt", os.path.basename(s_file))
        try:
            sample = tmp_s.group(1)
            kmer_all[sample] = hits_look_up(s_file)
        except AttributeError:
            continue

    # save results
    os.makedirs(out_folder, exist_ok=True)
    f_name = os.path.join(out_folder, "kmer_hits.json")
    save_data_to_json(kmer_all, f_name)

    print("job completed")
