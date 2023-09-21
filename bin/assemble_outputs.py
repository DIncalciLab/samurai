#!/usr/bin/env python3

import argparse

import janitor  # noqa:F401
import pandas as pd
from natsort import natsort_keygen


def parse_output(filename):

    # The rows past the third and up to the 17th contain the information
    # we require
    dataframe = pd.read_table(filename, skiprows=3,
                              nrows=14, sep=r":\t|:\s+",
                              engine="python")

    dataframe = dataframe.loc[["Tumor Fraction", "Ploidy",
                               "GC-Map correction MAD"]]

    return dataframe


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--samplesheet",
                        help="Sample sheet with metadata")
    parser.add_argument("source_files", nargs="+")
    options = parser.parse_args()

    sources = options.source_files

    combined_df = pd.concat([parse_output(item) for item in sources],
                            axis=1).transpose()
    combined_df.index.name = "samplename"

    if options.samplesheet is not None:

        samplesheet = (
                pd.read_csv(options.samplesheet, index_col=["samplename"])
                .select_columns(["patient", "type"])
        )
        # Avoid casts to other types, which will cause silent merge failures
        samplesheet.index = samplesheet.index.astype("str")
        combined_df = (
            combined_df
            .merge(samplesheet, left_index=True, right_index=True)
            .reset_index()
            .reorder_columns(["patient", "samplename", "type"])
            .rename_column("GC-Map correction MAD", "MAD")
            .clean_names()
            .change_type("tumor_fraction", float)
            .change_type("mad", float)
            .transform_column("tumor_fraction", lambda x: x.mul(100),
                              elementwise=False)
            .sort_values(by=["patient", "samplename"], key=natsort_keygen())
            .rename_column("patient", "patient_name")
            .reorder_columns(["samplename", "patient_name", "type", "ploidy",
                              "tumor_fraction", "mad"])
        )

    else:
        combined_df = combined_df.reset_index()

    combined_df.to_csv("ploidy_summary_mqc.txt", sep="\t", index=False)
    combined_df.to_excel("ichorCNA_summary.xlsx", sheet_name="Summary",
                         index=False)


if __name__ == "__main__":
    main()
