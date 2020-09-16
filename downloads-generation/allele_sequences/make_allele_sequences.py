"""
Generate allele sequences for pan-class II models.

Additional dependency: biopython
"""
from __future__ import print_function

import sys
import argparse

import pandas

import Bio.SeqIO  # pylint: disable=import-error

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "aligned_fasta",
    help="Aligned sequences")

parser.add_argument(
    "--reference-allele",
    required=True,
    help="Allele to use for position numbering")

parser.add_argument(
    "--out-csv",
    help="Result file")


def run():
    args = parser.parse_args(sys.argv[1:])
    print(args)

    allele_to_sequence = {}
    reader = Bio.SeqIO.parse(args.aligned_fasta, "fasta")
    for record in reader:
        name = record.description.split()[1]
        print(record.name, record.description)
        allele_to_sequence[name] = str(record.seq)

    allele_to_sequence = pandas.Series(allele_to_sequence).sort_index()
    print("Read %d aligned sequences" % len(allele_to_sequence))

    reference = allele_to_sequence[args.reference_allele]
    print("Using reference", args.reference_allele, reference)

    df = pandas.DataFrame(index=allele_to_sequence.index)

    current_number = 1
    for (i, reference_char) in enumerate(reference):
        if current_number not in df.columns:
            df[current_number] = ""

        df[current_number] += allele_to_sequence.str.get(i)
        if reference_char != '-':
            current_number += 1

    df = df.applymap(lambda s: s.replace("-", "X"))
    print(df)

    df.to_csv(args.out_csv, index=True)
    print("Wrote [%d alleles]: %s" % (len(df), args.out_csv))


if __name__ == '__main__':
    run()
