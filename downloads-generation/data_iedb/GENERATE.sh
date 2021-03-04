#!/bin/bash
#
# Download latest database snapshot from IEDB.
#
set -e
set -x

DOWNLOAD_NAME=data_iedb
SCRATCH_DIR=${TMPDIR-/tmp}/mhc2flurry-downloads-generation
SCRIPT_ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

mkdir -p "$SCRATCH_DIR"
rm -rf "$SCRATCH_DIR/$DOWNLOAD_NAME"
mkdir "$SCRATCH_DIR/$DOWNLOAD_NAME"

# Send stdout and stderr to a logfile included with the archive.
exec >  >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt")
exec 2> >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt" >&2)

# Log some environment info
date

cd $SCRATCH_DIR/$DOWNLOAD_NAME

wget -q https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip -O mhc_ligand_full.zip
wget -q http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip -O tcell_full_v3.zip

# We also download Uniprot viral sequences.
#
# UniprotKB search:
#    proteome:(excluded:no taxonomy:"Viruses [10239]") host:"Homo sapiens (Human) [9606]"
# link to human readable:
# https://www.uniprot.org/uniprot/?query=proteome%3A%28excluded%3Ano+taxonomy%3A%22Viruses+%5B10239%5D%22%29+AND+host%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score
wget -q 'https://www.uniprot.org/uniprot/?query=proteome:(excluded:no%20taxonomy:%22Viruses%20[10239]%22)%20AND%20host:%22Homo%20sapiens%20(Human)%20[9606]%22&format=fasta&sort=score' -O viruses.uniprot.fasta
gzip viruses.uniprot.fasta

unzip mhc_ligand_full.zip
rm mhc_ligand_full.zip
bzip2 mhc_ligand_full.csv

unzip tcell_full_v3.zip
rm tcell_full_v3.zip
bzip2 tcell_full_v3.csv

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
