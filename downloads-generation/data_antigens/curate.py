"""
Determine source proteins for antigens in IEDB.
"""
import sys
import os
import argparse

import pandas

import shellinford
import tqdm

from mhc2flurry.fasta import read_fasta_to_dataframe

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--data-iedb",
    metavar="tcell_full_v3.csv",
    help="Path to IEDB-style T cell epitope data")

parser.add_argument(
    "--data-sequences",
    metavar="SEQ.fasta",
    help="Path to IEDB-style T cell epitope data (tcell_full_v3.csv)")

parser.add_argument(
    "--out-csv",
    required=True,
    help="Combined result file")


def run():
    args = parser.parse_args(sys.argv[1:])

    epitopes_df = pandas.read_csv(args.data_iedb, skiprows=1)
    print("Read epitopes", *epitopes_df.shape)
    print(epitopes_df)

    sequences_df = read_fasta_to_dataframe(
        args.data_sequences, full_descriptions=True)
    print("Read sequences", len(sequences_df))
    print(sequences_df)

    # Prefer longest first
    sequences_df["length"] = sequences_df.sequence.str.len()
    sequences_df = sequences_df.sort_values("length", ascending=False)

    print("Building index")
    fm = shellinford.FMIndex()
    fm.build(sequences_df.sequence.tolist())
    print("Done")

    match_indices = {}
    peptides = epitopes_df.Description.dropna().str.upper().unique()
    for peptide in tqdm.tqdm(peptides):
        docs = sorted(fm.search(peptide), key=lambda doc: doc.doc_id)
        if len(docs) == 0:
            result = None
        else:
            result = docs[0].doc_id
        match_indices[peptide] = result

    result = pandas.DataFrame(index=peptides)
    result["antigen_index"] = result.index.map(match_indices)
    result = result.dropna().copy()
    result["antigen_name"] = result["antigen_index"].map(
        sequences_df.sequence_id)
    result["antigen_sequence"] = result["antigen_index"].map(
        sequences_df.sequence)

    print("Generated result", *result.shape)
    print(result)

    result.to_csv(args.out_csv, index=True)
    print("Wrote", args.out_csv)


if __name__ == '__main__':
    run()
