#!/bin/bash
#
# Class II structures from PDB
#
# Requires: curl
#
set -e
set -x

DOWNLOAD_NAME=data_pdb
SCRATCH_DIR=${TMPDIR-/tmp}/mhc2flurry-downloads-generation
SCRIPT_ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR=$(dirname "$SCRIPT_ABSOLUTE_PATH")
export PYTHONUNBUFFERED=1

mkdir -p "$SCRATCH_DIR"
rm -rf "$SCRATCH_DIR/$DOWNLOAD_NAME"
mkdir "$SCRATCH_DIR/$DOWNLOAD_NAME"

# Send stdout and stderr to a logfile included with the archive.
exec >  >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt")
exec 2> >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt" >&2)

# Log some environment info
date
pip freeze
git status
which clustalo
clustalo --version

cd $SCRATCH_DIR/$DOWNLOAD_NAME
cp $SCRIPT_DIR/batch_download.sh .
cp $SCRIPT_DIR/*.json .
cp $SCRIPT_ABSOLUTE_PATH .

curl https://search.rcsb.org/rcsbsearch/v1/query --data-urlencode "json@pdb_query.json" -G -o results.json

# TODO: make the query programmatically to include multiple alleles (not just DRA*01:01 as now).
# TODO: download the entries.

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
