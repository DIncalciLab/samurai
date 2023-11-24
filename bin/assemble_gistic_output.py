#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import pandas as pd
from gistic import GisticTable

def rearrange_output(gistic_file):

    folder_path = os.path.dirname(gistic_file)
    gistic_data = GisticTable.from_path(folder_path)
    
    return gistic_data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gistic_file", help="GISTIC output file from which get the directory.")
    options = parser.parse_args()

    gistic_file = options.gistic_file

    gistic_data = rearrange_output(gistic_file)
    all_lesions = gistic_data.all_lesions
    log_r_table = gistic_data.log2ratio_table
    cn_states = gistic_data.copy_number_state_table
    
    # FIX ME:
    # The gene table is 'None', but it exists if I run the function on same data with the container alone
    # all_genes_table = gistic_data.all_genes_table 

    all_lesions.to_csv("gistic_lesions.txt", sep="\t", index=False)
    # all_genes_table.to_csv("gistic_genes.txt", sep="\t", index=False)
    log_r_table.to_csv("gistic_log2R.txt", sep="\t", index=False)
    cn_states.to_csv("gistic_cn_states.txt", sep="\t", index=False)
    


if __name__ == "__main__":
    main()