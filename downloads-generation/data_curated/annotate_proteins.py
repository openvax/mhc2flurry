"""
Given a CSV where some column indicates peptides, add a column indicating which
protein(s) from some specified proteome contain that peptide.
"""

import argparse
import pandas
import time
import tqdm

import sys

import shellinford

from mhc2flurry.fasta import read_fasta_to_dataframe

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "input_csv",
    metavar="CSV",
    help="Input file")
parser.add_argument(
    "proteome",
    metavar="FASTA",
    help="Fasta proteome to search")
parser.add_argument(
    "--peptide-column",
    default="peptide",
    help="Name of column that gives peptides. Default: %(default)s")
parser.add_argument(
    "--protein-column",
    default="proteins",
    help="Name of column to write proteins. Default: %(default)s")
parser.add_argument(
    "--full-descriptions",
    default=False,
    action="store_true",
    help="Write the full protein descriptions, not just the IDs.")
parser.add_argument(
    "--join-character",
    default=" ",
    help="Separator to use between protein names. Default: '%(default)s'")
parser.add_argument(
    "--out-csv",
    required=True,
    metavar="CSV",
    help="Result file")


def run():
    args = parser.parse_args(sys.argv[1:])

    df = pandas.read_csv(args.input_csv)
    proteome_df = read_fasta_to_dataframe(
        args.proteome, full_descriptions=args.full_descriptions)

    print("Building FM index")
    start = time.time()
    fm = shellinford.FMIndex()
    fm.build(proteome_df.sequence.tolist())
    print("Built index of %d sequences in %0.3f sec." % (
        len(proteome_df), time.time() - start))

    result = []
    for peptide in tqdm.tqdm(df[args.peptide_column].values):
        matches = [item.doc_id for item in fm.search(peptide)]
        names = args.join_character.join(
            proteome_df.loc[matches, "sequence_id"].values)
        result.append(names)

    df[args.protein_column] = result
    df.to_csv(args.out_csv, index=False)
    print("Wrote", args.out_csv)


if __name__ == '__main__':
    run()
