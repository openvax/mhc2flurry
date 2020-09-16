"""
Filter and combine class II sequence fastas.
"""
from __future__ import print_function

import sys
import argparse

import mhcnames

import Bio.SeqIO  # pylint: disable=import-error


def normalize(s, disallowed=["MIC", "HFE"]):
    if any(item in s for item in disallowed):
        return None
    try:
        return mhcnames.normalize_allele_name(s, infer_class2_pair=False)
    except:
        while s:
            s = ":".join(s.split(":")[:-1])
            try:
                return mhcnames.normalize_allele_name(s, infer_class2_pair=False)
            except:
                pass

        print("Couldn't parse", s)
        return None


parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "fastas",
    nargs="+",
    help="Unaligned fastas")

parser.add_argument(
    "--kind",
    required=True,
    choices=("alpha", "beta"),
    help="Chain")

parser.add_argument(
    "--out",
    required=True,
    help="Fasta output")

min_lengths = {
    "alpha": 200,
    "beta": 200,
}


def run():
    args = parser.parse_args(sys.argv[1:])
    print(args)

    min_length = min_lengths[args.kind]

    output_records = []
    seen = set()
    sequences = set()

    input_records = []
    for fasta in args.fastas:
        reader = Bio.SeqIO.parse(fasta, "fasta")
        input_records.extend(reader)

    # Iterate longest records first so that when multiple records have the
    # same two digit normalized allele, we use the longest one.
    for record in sorted(input_records, key=lambda r: len(r.seq), reverse=True):
        name = record.description.split()[1]
        name = normalize(name)
        if name in seen:
            continue
        if len(record.seq) < min_length:
            print("Skipping due to short length", name, record.description)
            continue
        seen.add(name)
        sequences.add(record.seq)
        record.description = name + " " + record.description
        output_records.append(record)

    with open(args.out, "w") as fd:
        Bio.SeqIO.write(output_records, fd, "fasta")

    print("Wrote %d / %d [%d unique] sequences: %s" % (
        len(output_records), len(input_records), len(sequences), args.out))


if __name__ == '__main__':
    run()
