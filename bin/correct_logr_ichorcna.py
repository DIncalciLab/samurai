#!/usr/bin/env python

from pathlib import Path

import polars as pl
import typer

from typing import Optional
from typing_extensions import Annotated

__version__ = "0.1.0"


def version_callback(value: bool):
    if value:
        print(f"correct_logr_ichorcna {__version__}")
        raise typer.Exit()


def main(
    seg: Path = typer.Option(
        default=...,
        help="Segmentation file from which extract information",
    ),
    ploidy: Path = typer.Option(
        default=...,
        help="Ploidy summary file from which extract information",
    ),
    version: bool | None = typer.Option(
        default=None, callback=version_callback
    ),
):

    segments = pl.scan_csv(seg, separator="\t").rename({"ID": "sample"})
    ploidy_df = pl.scan_csv(ploidy, separator="\t")

    result = (
        segments.join(
            ploidy_df.select(pl.col("samplename").alias("sample"), "Ploidy"),
            on="sample",
            how="left",
        )
        .with_columns(
            pl.col("logR_Copy_Number")
            .truediv(pl.col("Ploidy"))
            .clip(lower_bound=2e-8)
            .log(base=2)
            .clip(lower_bound=-0.5)
            .alias("adj.seg")
        )
        .select("sample", "chromosome", "start", "end", "num.mark", "adj.seg")
    )

    result.collect().write_csv(
        "segments_logR_corrected_gistic.seg",
        separator="\t",
        quote_style="necessary",
    )


if __name__ == "__main__":
    typer.run(main)
