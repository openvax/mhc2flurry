#!/bin/bash
#
# Class II allele sequences
#
# Requires: clustalo, wget
#
set -e
set -x

DOWNLOAD_NAME=allele_sequences
SCRATCH_DIR=${TMPDIR-/tmp}/mhcflurryii-downloads-generation
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
cp $SCRIPT_DIR/make_allele_sequences.py .
cp $SCRIPT_DIR/filter_sequences.py .
cp $SCRIPT_ABSOLUTE_PATH .

# Human

# Alpha chain
mkdir alpha
cd alpha
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPA1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPA2_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DQA1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DQA2_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DRA_prot.fasta
cd ..

# Beta chain
mkdir beta
cd beta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPB1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPB2_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DQB1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DRB_prot.fasta
cd ..

python filter_sequences.py alpha/*.fasta --kind alpha --out alpha.fasta
python filter_sequences.py beta/*.fasta --kind beta --out beta.fasta

time clustalo -i alpha.fasta -o alpha.aligned.fasta
time clustalo -i beta.fasta -o beta.aligned.fasta

time python make_allele_sequences.py \
    alpha.aligned.fasta \
    --reference-allele HLA-DRA1*01:01 \
    --out-csv alpha.csv

time python make_allele_sequences.py \
    beta.aligned.fasta \
    --reference-allele HLA-DRB1*01:01 \
    --out-csv beta.csv

# Cleanup
gzip -f alpha.fasta
gzip -f alpha.aligned.fasta
gzip -f beta.fasta
gzip -f beta.aligned.fasta

for i in $(ls "*/*.fasta")
do
    gzip -f $i
done

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
