"""
Determine source proteins for antigens in IEDB.
"""
import sys
import os
import json
import time
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
    "--max-epitopes",
    metavar="N",
    type=int,
    help="Process first N epitopes (for debugging)")

parser.add_argument(
    "--max-sequences",
    metavar="N",
    type=int,
    help="Subsample to N sequences (for debugging)")

parser.add_argument(
    "--min-sequence-length",
    metavar="N",
    type=int,
    default=20,
    help="Discard sequences shorter than N. Default: %(default)s")

parser.add_argument(
    "--out-csv",
    required=True,
    help="Combined result file")


def run():
    args = parser.parse_args(sys.argv[1:])

    epitopes_df = pandas.read_csv(
        args.data_iedb, skiprows=1, nrows=args.max_epitopes)
    print("Read epitopes", *epitopes_df.shape)
    print(epitopes_df)

    sequences_df = read_fasta_to_dataframe(
        args.data_sequences, full_descriptions=True)
    if args.min_sequence_length:
        sequences_df = sequences_df.loc[
            sequences_df.sequence.str.len() >= args.min_sequence_length
        ]
    if args.max_sequences:
        sequences_df = sequences_df.sample(n=args.max_sequences)
    print("Read sequences", len(sequences_df))
    print(sequences_df)

    # Prefer longest first
    sequences_df["length"] = sequences_df.sequence.str.len()
    sequences_df = sequences_df.sort_values(
        "length", ascending=False).drop_duplicates(
        "sequence").reset_index(drop=True)

    print("Building index")
    start = time.time()
    fm = shellinford.FMIndex()
    fm.build(sequences_df.sequence.tolist())
    print("Done [%0.2f sec.]" % (time.time() - start))

    match_indices = {}
    peptides = epitopes_df.Description.dropna().str.upper().unique()
    for peptide in tqdm.tqdm(peptides):
        docs = fm.search(peptide)
        match_indices[peptide] = sorted([d.doc_id for d in docs])

    result = pandas.DataFrame(index=peptides)
    result.index.name = "peptide"
    result["antigen_indices"] = result.index.map(match_indices)
    result = result.loc[
        result.antigen_indices.str.len() > 0
    ].copy()
    result["antigen_index"] = result["antigen_indices"].str.get(0)

    result["antigen_name"] = result["antigen_index"].map(
        sequences_df.sequence_id)
    result["antigen_sequence"] = result["antigen_index"].map(
        sequences_df.sequence)

    result["other_antigen_names"] = result["antigen_indices"].map(
        lambda some: sequences_df.loc[some[1:], "sequence_id"].tolist())
    result["other_antigen_names"] = result["other_antigen_names"].map(json.dumps)

    del result['antigen_index']
    del result['antigen_indices']

    print("Generated result", *result.shape)
    print(result)

    result.to_csv(args.out_csv, index=True)
    print("Wrote", args.out_csv)


if __name__ == '__main__':
    run()
