# From a PDB results json, print out a comma separated list of PDB IDs

import argparse
import sys
import json

parser = argparse.ArgumentParser()
parser.add_argument("results", metavar="JSON")
parser.add_argument("out", metavar="FILE")
args = parser.parse_args(sys.argv[1:])

parsed = json.load(open(args.results))
print("Loaded %d results" % len(parsed['result_set']))

print("First result")
print(parsed['result_set'][0])

print("Last result")
print(parsed['result_set'][-1])

with open(args.out, "w") as fd:
    identifiers = [entry['identifier'] for entry in parsed['result_set']]
    fd.write(",".join(identifiers))
    fd.write("\n")
print("Wrote: ", args.out)
