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
    nargs="+",
    help="Aligned sequences")

parser.add_argument(
    "--reference-allele",
    help="Allele to use for position numbering. If not specified, a single "
    "column is written out with the full sequence. If specified, multiple "
    "columns are written giving the sequence at each position in the reference "
    "allele.")

parser.add_argument(
    "--out-csv",
    help="Result file")


def run():
    args = parser.parse_args(sys.argv[1:])
    print(args)

    allele_to_sequence = {}
    for fasta in args.aligned_fasta:
        reader = Bio.SeqIO.parse(fasta, "fasta")
        for record in reader:
            name = record.description.split()[1]
            print(record.name, record.description)
            allele_to_sequence[name] = str(record.seq)

    allele_to_sequence = pandas.Series(allele_to_sequence).sort_index()
    print("Read %d aligned sequences" % len(allele_to_sequence))

    if args.reference_allele:
        # Output multiple columns for each position in the reference allele
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
    else:
        # Output just a single column with the whole sequence
        df = allele_to_sequence.to_frame()
        df.index.name = 'allele'
        df.columns = ['sequence']

    df = df.applymap(lambda s: s.replace("-", "X"))
    print(df)

    df.to_csv(args.out_csv, index=True)
    print("Wrote [%d alleles]: %s" % (len(df), args.out_csv))


if __name__ == '__main__':
    run()
