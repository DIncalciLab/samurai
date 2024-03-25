#!/usr/bin/env python3

import argparse
from typing import Optional

import janitor  # noqa: F401
import pandas as pd
from pybedtools import BedTool


def add_sample_data(dataframe: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    dataframe = (
        dataframe.add_column("id", sample_id)
        .reorder_columns(["id"])
        .reorder_columns(
            ["id", "chrom", "start", "end", "num.mark", "seg.mean.adj", "call"]
        )
    )

    return dataframe


def convert_to_seg(
    segments: pd.DataFrame,
    bins: pd.DataFrame,
    aberrations: pd.DataFrame,
    sample: Optional[str] = None,
) -> pd.DataFrame:
    segments_bed = BedTool.from_dataframe(segments)
    bins_bed = BedTool.from_dataframe(bins)
    aberrations_bed = BedTool.from_dataframe(aberrations)

    result = segments_bed.intersect(bins_bed, F=1, c=True).to_dataframe(
        disable_auto_names=True,
        names=["chrom", "start", "end", "seg.mean", "zscore", "num.mark"],
    )

    if aberrations.empty:
        # No alterations, return segments with zeroed log2ratios
        result = result.transform_column("seg.mean", lambda x: 0)
        result = (
            result.reorder_columns(
                ["chrom", "start", "end", "num.mark", "seg.mean", "zscore"]
            )
            .rename_column("seg.mean", "seg.mean.adj")
            .remove_columns("zscore")
            .add_column("call", "neut")
        )

        if sample is not None:
            result = add_sample_data(result, sample)

        return result

    result = BedTool.from_dataframe(result.remove_columns("zscore")).intersect(
        aberrations_bed, loj=True
    )

    result = (
        result.to_dataframe(
            disable_auto_names=True,
            names=[
                "chrom",
                "start",
                "end",
                "seg.mean",
                "num.mark",
                "a_chr",
                "a_start",
                "a_end",
                "a_logr",
                "a_zscore",
                "call",
            ],
        )
        .process_text(
            "call", string_function="replace", pat=".", repl="neut", regex=False
        )
        .select_columns(["chrom", "start", "end", "num.mark", "seg.mean", "call"])
        #  Zero logR if the call is neutral
        .join_apply(
            lambda x: 0 if x.call == "neut" else x["seg.mean"],
            new_column_name="seg.mean.adj",
        )
        .remove_columns(["seg.mean"])
    )

    if sample is not None:
        result = add_sample_data(result, sample)

    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--id", metavar="SAMPLEID", help="Add SAMPLEID to the converted file"
    )
    parser.add_argument("segments", help="File containing segments")
    parser.add_argument("bins", help="File containing bins")
    parser.add_argument("aberrations", help="File containing aberration calls")
    parser.add_argument("destination", help="File to save results to")
    parser.add_argument("destination_gistic", help="File to save results to")

    options = parser.parse_args()
    segments = pd.read_table(options.segments)
    bins = pd.read_table(options.bins)
    aberrations = pd.read_table(options.aberrations)
    destination = options.destination
    destination_gistic = options.destination_gistic

    converted_df = convert_to_seg(segments, bins, aberrations, sample=options.id)

    converted_df.to_csv(destination, sep="\t", index=False)

    # Remove the 'call' column
    gistic_df = converted_df.remove_columns("call").rename_column("seg.mean.adj", "seg.mean")

    gistic_df.to_csv(destination_gistic, sep="\t", index=False)


if __name__ == "__main__":
    main()
