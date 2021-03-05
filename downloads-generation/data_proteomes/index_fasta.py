"""
Write a shellinford index for a fasta.
"""

import argparse
import time
import sys

import shellinford

from mhc2flurry.fasta import read_fasta_to_dataframe

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "input",
    metavar="FASTA",
    help="Input file")
parser.add_argument(
    "output",
    metavar="FM",
    help="Output file")


def run():
    args = parser.parse_args(sys.argv[1:])

    df = read_fasta_to_dataframe(args.input)
    print("Read")
    print(df)

    print("Building FM index")
    start = time.time()
    fm = shellinford.FMIndex()
    fm.build(df.sequence.tolist())
    print("Built index of %d sequences in %0.3f sec." % (
        len(df), time.time() - start))

    print("Writing index")
    fm.write(args.output)
    print("Wrote", args.output)



if __name__ == '__main__':
    run()
