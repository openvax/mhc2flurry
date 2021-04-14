"""
Curate IEDB T cell epitopes. Currently this doesn't do much except rename the
peptide column from "Description" to "peptide".
"""
import sys
import argparse

import pandas

from mhc2flurry.amino_acid import COMMON_AMINO_ACIDS

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--data-iedb",
    metavar="tcell_full_v3.csv",
    help="Path to IEDB-style T cell epitope data")

parser.add_argument(
    "--max-epitopes",
    metavar="N",
    type=int,
    help="Process first N epitopes (for debugging)")

parser.add_argument(
    "--out-csv",
    required=True,
    help="Result file")


def run():
    args = parser.parse_args(sys.argv[1:])

    epitopes_df = pandas.read_csv(
        args.data_iedb, skiprows=1, nrows=args.max_epitopes)
    print("Read epitopes", *epitopes_df.shape)
    print(epitopes_df)

    epitopes_df.insert(0, "peptide", epitopes_df.Description)
    aa_regex = "^[%s]+$" % "".join(sorted(COMMON_AMINO_ACIDS))

    epitopes_df = epitopes_df.loc[
        epitopes_df.peptide.str.match(aa_regex) &
        (epitopes_df.peptide.str.len() >= 5)
    ]

    print("Epitopes with valid peptides", len(epitopes_df))

    print("Generated result", *epitopes_df.shape)
    print(epitopes_df)

    epitopes_df.to_csv(args.out_csv, index=False)
    print("Wrote", args.out_csv)


if __name__ == '__main__':
    run()
