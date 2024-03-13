#!/usr/bin/env python3

import argparse

import janitor  # noqa: F401
from janitor.io import read_csvs
import pandas as pd


def extract_data(filename: str, source: pd.DataFrame) -> pd.DataFrame:
    samplename = filename.replace("_aberrations.bed", "")

    df = source.add_column("id", samplename).reorder_columns(["id"])

    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_files", nargs="+")

    options = parser.parse_args()

    source_data = read_csvs(options.source_files, sep="\t", separate_df=True)

    final_df = pd.concat(
        [extract_data(str(filename), df) for filename, df in source_data.items()]
    )

    final_df.to_csv("all_aberrations.txt", index=False, sep="\t")


if __name__ == "__main__":
    main()
