#!/bin/bash
#
# Download viral, human, and mouse proteomes from Uniprot
#
set -e
set -x

DOWNLOAD_NAME=data_proteomes
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

###### Human protein sequences
# For now we are using the "one protein sequence per gene" version provided by
# Uniprot. We'll see if we end up missing a lot of peptides this way (i.e. by
# excluding alternate isoforms).
wget -q 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz' -O human.uniprot.one_per_gene.fasta.gz

# We also download ID mapping to Ensembl and other databases:
wget -q 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz' -O human.uniprot_id_mapping.dat.gz

###### Mouse protein sequences
# Also using the "one proein sequence per gene" version.
wget -q 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz' -O mouse.uniprot.one_per_gene.fasta.gz
wget -q 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz' -O mouse.uniprot_id_mapping.dat.gz

###### Viral protein sequences
# UniprotKB search:
#    proteome:(excluded:no taxonomy:"Viruses [10239]") host:"Homo sapiens (Human) [9606]"
# link to human readable:
# https://www.uniprot.org/uniprot/?query=proteome%3A%28excluded%3Ano+taxonomy%3A%22Viruses+%5B10239%5D%22%29+AND+host%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score
wget -q 'https://www.uniprot.org/uniprot/?query=proteome:(excluded:no%20taxonomy:%22Viruses%20[10239]%22)%20AND%20host:%22Homo%20sapiens%20(Human)%20[9606]%22&format=fasta&sort=score' -O viral.uniprot.fasta
gzip viral.uniprot.fasta

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
