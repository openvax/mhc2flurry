# Given a set of PDB .cif.gz files, write out a fasta with the sequences of
# each chain. This will be used to align MHC II PDB structures against
# sequences from IMDB and other sources.

import argparse
import sys
import json
import os
import glob

import atomium

parser = argparse.ArgumentParser()
parser.add_argument(
    "input", metavar="JSON", help='Director of .cif.gz files')
parser.add_argument("out", metavar="FILE.fasta", help="Out fasta file")
args = parser.parse_args(sys.argv[1:])

print(args)

files = glob.glob(args.input + "/*.cif.gz")
print("Found %d files" % len(files))

with open(args.out, "w") as fd:
    for file in files:
        structure = atomium.open(file)
        for chain in structure.model.chains():
            fd.write(">%s_%s %s\n" % (
                structure.code, chain.id, os.path.basename(file)))
            fd.write("".join(c.code for c in chain.residues()))
            fd.write("\n")
print("Wrote: ", args.out)
