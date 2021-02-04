#!/bin/bash
#
# Full protein sequences for epitopes in IEDB
# Requirement: shellinford
set -e
set -x

DOWNLOAD_NAME=data_antigens
SCRATCH_DIR=${TMPDIR-/tmp}/mhc2flurry-downloads-generation
SCRIPT_ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR=$(dirname "$SCRIPT_ABSOLUTE_PATH")

mkdir -p "$SCRATCH_DIR"
rm -rf "$SCRATCH_DIR/$DOWNLOAD_NAME"
mkdir "$SCRATCH_DIR/$DOWNLOAD_NAME"

# Send stdout and stderr to a logfile included with the archive.
exec >  >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt")
exec 2> >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt" >&2)

# Log some environment info
date

cd $SCRATCH_DIR/$DOWNLOAD_NAME
cp $SCRIPT_DIR/curate.py .

time python curate.py \
    --data-iedb \
        "$(mhc2flurry-downloads path data_iedb)/tcell_full_v3.csv.bz2" \
    --data-sequences \
        "$(mhc2flurry-downloads path data_iedb)/viruses.uniprot.fasta.gz" \
    --out-csv epitope_to_sequence.csv

bzip2 epitope_to_sequence.csv

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
