#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Format ichorCNA .seg.txt file ")
    parser.add_argument("-s", "--segmentation_file", required=True, help="Input ichorCNA .seg file")
    parser.add_argument("-o", "--output", default=None, help="Output formatted .seg file")
    args = parser.parse_args()

    # Create file
    df = pd.read_csv(args.segmentation_file, sep="\t", dtype=str)

    # Check columns in .seg.txt file
    required_cols = ["chrom", "start", "end", "logR_Copy_Number", "ID"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"There mandatory columns are missing {args.segmentation_file}: {missing}")

    # Rename columns as required for CX signatures
    df_out = (
        df.loc[:, required_cols]
          .rename(columns={
              "chrom": "chromosome",
              "logR_Copy_Number": "segVal",
              "ID": "sample"
          })
    )

    # Convert to numeric 
    for col in ["start", "end", "segVal"]:
        df_out[col] = pd.to_numeric(df_out[col], errors="coerce")

    # Remove NAs
    df_out = df_out.dropna(subset=["segVal"])

    # Output file
    out_file = args.output 

    # Save formatted segmentation file
    df_out.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
