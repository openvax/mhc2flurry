"""
Given a CSV where some column indicates peptides, add a column indicating which
protein(s) from some specified proteome contain that peptide.
"""

import argparse
import time
import sys

import tqdm
import pandas
import numpy
import shellinford

from mhc2flurry.fasta import read_fasta_to_dataframe

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "reference",
    metavar="FASTA",
    help="Fasta proteome to search.")
parser.add_argument(
    "--annotate",
    action="append",
    default=[],
    nargs=2,
    metavar="CSV",
    help="Input and output file pairs. Specify this argument multiple times "
    "to process multiple input files, each of which will be written to its "
    "respective output file. The output file can be specified as '-' to "
    "overwrite the input file.")
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
    "--fm-index-suffix",
    metavar="SUFFIX",
    help="Use a pre-existing fm index found by concatenating SUFFIX onto each "
    "input fasta filename.")


def run():
    args = parser.parse_args(sys.argv[1:])

    peptides = set()
    input_filename_df_and_output_filename = []

    for (input, output) in args.annotate:
        if output.strip() == "-":
            output = input
        df = pandas.read_csv(input)
        print("Read peptides", input)
        print(df)
        input_filename_df_and_output_filename.append((input, df, output))
        peptides.update(df[args.peptide_column].unique())

    print("Read %d peptides to annotate" % len(peptides))

    proteome_df = read_fasta_to_dataframe(
        args.reference, full_descriptions=args.full_descriptions)

    print("Read proteome:")
    print(proteome_df)

    fm = shellinford.FMIndex()
    start = time.time()
    if args.fm_index_suffix:
        name = args.reference + args.fm_index_suffix
        print("Using pre-existing fm index", name)
        fm.read(name)
        print("Read in %0.3f sec." % (time.time() - start))
    else:
        print("Building FM index")
        fm.build(proteome_df.sequence.tolist())
        print("Built index of %d sequences in %0.3f sec." % (
            len(proteome_df), time.time() - start))

    print("Annotating peptides")
    peptide_to_matches = {}
    for peptide in tqdm.tqdm(peptides):
        matches = [item.doc_id for item in fm.search(peptide)]
        names = args.join_character.join(
            proteome_df.loc[matches, "sequence_id"].values)
        peptide_to_matches[peptide] = names

    print("Writing files")
    for (input, df, output) in input_filename_df_and_output_filename:
        print(input)
        df[args.protein_column] = df[args.peptide_column].map(
            peptide_to_matches)
        df.to_csv(output, index=False)
        print("Wrote", output)


if __name__ == '__main__':
    run()
