# Assign PDB sequences (searched by mmseqs against IMGT sequences)
# to alpha vs beta based on mmseqs results

import argparse
import sys
import pandas
import os

from mhc2flurry.fasta import read_fasta_to_dataframe

parser = argparse.ArgumentParser()
parser.add_argument(
    "pdb_sequences",
    metavar="FASTA",
    help='PDB sequences')
parser.add_argument(
    "search_results",
    metavar="TXT",
    help='mmseqs search results')
parser.add_argument(
    "--mmseqs-output-format",
    metavar="A,B,C",
    required=True,
    help='mmseqs output format (comma separated list of fields)')
parser.add_argument(
    "--out-alpha",
    metavar="FASTA",
    help='Output file')
parser.add_argument(
    "--out-beta",
    metavar="FASTA",
    help='Output file')
args = parser.parse_args(sys.argv[1:])

print(args)

sequences_df = read_fasta_to_dataframe(args.pdb_sequences).set_index("sequence_id")

search_df = pandas.read_csv(
    args.search_results,
    names=args.mmseqs_output_format.split(","),
    sep=None)
search_df["kind"] = search_df.target.str.split(".").str.get(0)

df = search_df.loc[
    (search_df.qcov > 0.7) &
    (search_df.tcov > 0.5)
].sort_values("evalue").drop_duplicates("query").set_index("query")

print(df)
print("Breakdown by kind [should be equal or nearly equal]")
print(df.kind.value_counts())


def write_fasta(filename, sub_df):
    with open(filename, "w") as fd:
        for name, row in sub_df.iterrows():
            seq = sequences_df.loc[name].sequence
            fd.write(">pdb.%s\n" % name)
            fd.write(seq)
            fd.write("\n")
    print("Wrote", filename, "with", len(sub_df), "sequences")

if args.out_alpha:
    write_fasta(args.out_alpha, df.loc[df.kind == "alpha"])

if args.out_beta:
    write_fasta(args.out_beta, df.loc[df.kind == "beta"])
