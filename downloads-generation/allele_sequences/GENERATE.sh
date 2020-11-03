#!/bin/bash
#
# Class II allele sequences
#
# Requires: clustalo, wget
#
set -e
set -x

DOWNLOAD_NAME=allele_sequences
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
cp $SCRIPT_DIR/make_allele_sequences.py .
cp $SCRIPT_DIR/filter_sequences.py .
cp $SCRIPT_ABSOLUTE_PATH .

######## Human
# Human Alpha chain
mkdir -p alpha
cd alpha
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPA1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPA2_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DQA1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DQA2_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DRA_prot.fasta
cd ..

# Human Beta chain
mkdir -p beta
cd beta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPB1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DPB2_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DQB1_prot.fasta
wget -q ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/DRB_prot.fasta
cd ..

######## MOUSE
#
# See http://www.imgt.org/IMGTrepertoireMH/Proteins/tables/index.php?species=mouse&gene=MH2-AB
#
#
# Mouse Alpha chain: TODO
#
# Mouse Beta chain: H-2-AB alleles
# Commented out for now until we have A-2-AA sequences
#mkdir -p beta
#cd beta
#wget -q https://www.uniprot.org/uniprot/P14483.fasta  # MH2-AB*01
#wget -q https://www.uniprot.org/uniprot/P01921.fasta  # MH2-AB*03
#wget -q https://www.uniprot.org/uniprot/Q31135.fasta  # MH2-AB*06
#wget -q https://www.uniprot.org/uniprot/O78197.fasta  # MH2-AB*10
#wget -q https://www.uniprot.org/uniprot/O62867.fasta  # MH2-AB*12
#wget -q https://www.uniprot.org/uniprot/P06342.fasta  # MH2-AB*13
#wget -q https://www.uniprot.org/uniprot/O78197.fasta  # MH2-AB*10
#wget -q https://www.uniprot.org/uniprot/O62867.fasta  # MH2-AB*12
#wget -q https://www.uniprot.org/uniprot/P06342.fasta  # MH2-AB*13
#wget -q https://www.uniprot.org/uniprot/P06343.fasta  # MH2-AB*14
#wget -q https://www.uniprot.org/uniprot/Q31184.fasta  # MH2-AB*15
#wget -q https://www.uniprot.org/uniprot/Q31131.fasta  # MH2-AB*16
#wget -q https://www.uniprot.org/uniprot/P06345.fasta  # MH2-AB*17
#wget -q https://www.uniprot.org/uniprot/P06346.fasta  # MH2-AB*20
#wget -q https://www.uniprot.org/uniprot/O19470.fasta  # MH2-AB*22


python filter_sequences.py alpha/*.fasta --kind alpha --out alpha.fasta
python filter_sequences.py beta/*.fasta --kind beta --out beta.fasta

time clustalo -i "$(pwd)/alpha.fasta" -o "$(pwd)/alpha.aligned.fasta"
time clustalo -i "$(pwd)/beta.fasta" -o "$(pwd)/beta.aligned.fasta"

time python make_allele_sequences.py \
    "$(pwd)/alpha.aligned.fasta" \
    --reference-allele HLA-DRA*01:01 \
    --out-csv "$(pwd)/alpha.csv"

time python make_allele_sequences.py \
    "$(pwd)/beta.aligned.fasta" \
    --reference-allele HLA-DRB1*01:01 \
    --out-csv "$(pwd)/beta.csv"

time python make_allele_sequences.py \
    "$(pwd)/alpha.aligned.fasta" \
    "$(pwd)/beta.aligned.fasta" \
    --out-csv "$(pwd)/all.csv"

# Cleanup
gzip -f alpha.fasta
gzip -f alpha.aligned.fasta
gzip -f beta.fasta
gzip -f beta.aligned.fasta

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
