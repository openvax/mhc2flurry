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
cp $SCRIPT_DIR/*.py .
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

python filter_sequences.py alpha/*.fasta --kind alpha --out alpha.database.fasta
python filter_sequences.py beta/*.fasta --kind beta --out beta.database.fasta

# Generate PDB sequences
time python extract_pdb_sequences.py \
    "$(mhc2flurry-downloads path data_pdb)/structures" \
    "$(pwd)/pdb_sequences.fasta"

# Search PDB sequences against downloaded IMDB
cat alpha.database.fasta beta.database.fasta > alpha_and_beta.database.fasta
MMSEQS_OUTPUT_FORMAT="query,target,qaln,taln,qseq,tseq,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,evalue"
mmseqs easy-search \
    "$(pwd)/pdb_sequences.fasta" \
    "$(pwd)/alpha_and_beta.database.fasta" \
    "$(pwd)/pdb_search.m8" \
    "${TMPDIR-/tmp}" \
    --format-output "$MMSEQS_OUTPUT_FORMAT" \
    -s 1.0
time python assign_pdb_sequences_to_alpha_or_beta.py \
    "$(pwd)/pdb_sequences.fasta" \
    "$(pwd)/pdb_search.m8" \
    --mmseqs-output-format "$MMSEQS_OUTPUT_FORMAT" \
    --out-alpha "$(pwd)/alpha.pdb.fasta" \
    --out-beta "$(pwd)/beta.pdb.fasta"

cat alpha.pdb.fasta alpha.database.fasta > alpha.combined.fasta
cat beta.pdb.fasta beta.database.fasta > beta.combined.fasta

# Run clustalo to generate multiple sequence alignments
time clustalo -i "$(pwd)/alpha.combined.fasta" -o "$(pwd)/alpha.aligned.fasta" \
    --clustering-out cluster.alpha.aux
time clustalo -i "$(pwd)/beta.combined.fasta" -o "$(pwd)/beta.aligned.fasta" \
    --clustering-out cluster.beta.aux

time python make_pseudosequences.py \
    "$(pwd)/alpha.aligned.fasta" \
    "$(pwd)/beta.aligned.fasta" \
    "$(mhc2flurry-downloads path data_pdb)/structures" \
    --criteria 0.2 6.0 0.3 \
    --criteria 0.1 8.0 0.1 \
    --criteria 0.05 10.0 0.05 \
    --reference-allele HLA-DRA*01:01 HLA-DRB1*01:01 \
    --reference-structure 3QXD \
    --reference-structure 5KSU \
    --out-csv "$(pwd)/allele_sequences.csv" \
    --out-aux-dir "$(pwd)/aux-info"  # Extra info

# Cleanup
rm pdb_search.m8
for i in $(ls *.fasta)
do
    gzip -f "$i"
done

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
