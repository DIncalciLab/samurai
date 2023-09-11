#!/usr/bin/env python3

import argparse
from typing import Tuple, Hashable

from janitor.io import read_csvs  # noqa:F401
import pandas as pd

MAPPINGS = {
    "Copy number profile abnormality (CPA) score (doi: 10.1186/s13073-020-00735-4)": "cpa",
    "Median segment variance per bin (doi: 10.1093/nar/gky1263)": "msv",
    "Number of reads": "reads"
}


def extract_data(record) -> pd.DataFrame:

    filename, source = record
    source = (
        source
        .rename(index=MAPPINGS)
        .squeeze()
    )

    samplename = filename.replace("_statistics.txt", "")

    df = pd.DataFrame(index=[samplename], columns=["reads", "cpa", "msv"])
    df.index.name = "sample"

    df.loc[samplename, "cpa"] = source["cpa"]
    df.loc[samplename, "msv"] = source["msv"]
    df.loc[samplename, "reads"] = source["reads"]

    return df


def main() -> None:

    parser = argparse.ArgumentParser()
    parser.add_argument("source_files", nargs="+")

    options = parser.parse_args()

    # Parse the horrid format of WisecondorX

    source_data = read_csvs(
        options.source_files,
        header=None,
        skiprows=24,
        sep=r"(?<!doi): ",
        engine="python",
        index_col=0,
        separate_df=True,
        names=["data"],
    )

    final_df = pd.concat([extract_data(record)
                          for record in source_data.items()])

    final_df = (
        final_df
        .reset_index()
        .sort_values(by="cpa", ascending=False)
    )

    final_df.to_excel("summary.xlsx", index=False)


if __name__ == "__main__":
    main()