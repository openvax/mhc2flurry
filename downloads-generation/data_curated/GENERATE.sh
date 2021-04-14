#!/bin/bash
#
# Create "curated" training data, which combines an IEDB download with additional
# published data, removes unusable entries, normalizes allele name, and performs
# other filtering and standardization.
#
set -e
set -x

DOWNLOAD_NAME=data_curated
SCRATCH_DIR=${TMPDIR-/tmp}/mhc2flurry-downloads-generation
SCRIPT_ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR=$(dirname "$SCRIPT_ABSOLUTE_PATH")
export PYTHONUNBUFFERED=1

# Check that needed downloads are available.
mhc2flurry-downloads path data_published
mhc2flurry-downloads path data_iedb
mhc2flurry-downloads path data_proteomes

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
mhc2flurry-downloads info


cd $SCRATCH_DIR/$DOWNLOAD_NAME

cp $SCRIPT_DIR/curate_mhc_ligands.py .
cp $SCRIPT_DIR/curate_t_cell_epitopes.py .
cp $SCRIPT_DIR/curate_ms_by_pmid.py .
cp $SCRIPT_DIR/annotate_proteins.py .

# Curate IEDB T cell epitopes
time python curate_t_cell_epitopes.py \
    --data-iedb "$(mhc2flurry-downloads path data_iedb)/tcell_full_v3.csv.bz2" \
    --out-csv "$(pwd)/t_cell_epitopes.csv"


MS_DIR="$(mhc2flurry-downloads path data_published)/ms"
cp -r "$MS_DIR" .

EXPRESSION_DIR="$(mhc2flurry-downloads path data_published)/expression"
cp -r "$EXPRESSION_DIR" .

CURATE_BY_PMID_ARGS=""
for pmid in $(ls ms)
do
    CURATE_BY_PMID_ARGS+=$(echo --ms-item $pmid ms/$pmid/* ' ')
done
for item in $(ls expression)
do
    CURATE_BY_PMID_ARGS+=$(echo --expression-item $item expression/$item/* ' ')
done

time python curate_ms_by_pmid.py $CURATE_BY_PMID_ARGS \
    --ms-out ms.by_pmid.csv \
    --expression-out rna_expression.csv \
    --expression-metadata-out rna_expression.metadata.csv

rm -rf ms

time python curate_mhc_ligands.py \
    --data-iedb \
        "$(mhc2flurry-downloads path data_iedb)/mhc_ligand_full.csv.bz2" \
    --data-additional-ms "$(pwd)/ms.by_pmid.csv" \
    --out-csv curated_training_data.csv \
    --out-affinity-csv curated_training_data.affinity.csv \
    --out-mass-spec-csv curated_training_data.mass_spec.csv

time python curate_mhc_ligands.py \
    --data-iedb \
        "$(mhc2flurry-downloads path data_iedb)/mhc_ligand_full.csv.bz2" \
    --out-csv curated_training_data.no_additional_ms.csv

# Annotate human proteins
time python annotate_proteins.py \
    "$(mhc2flurry-downloads path data_proteomes)/human.uniprot.isoforms.fasta.gz" \
    --annotate "$(pwd)/ms.by_pmid.csv" - \
    --annotate "$(pwd)/curated_training_data.csv" - \
    --annotate "$(pwd)/curated_training_data.mass_spec.csv" - \
    --annotate "$(pwd)/curated_training_data.no_additional_ms.csv" - \
    --annotate "$(pwd)/t_cell_epitopes.csv" - \
    --fm-index-suffix .fm \
    --protein-column proteins_human

# Annotate mouse proteins
time python annotate_proteins.py \
    "$(mhc2flurry-downloads path data_proteomes)/mouse.uniprot.isoforms.fasta.gz" \
    --annotate "$(pwd)/ms.by_pmid.csv" - \
    --annotate "$(pwd)/curated_training_data.csv" - \
    --annotate "$(pwd)/curated_training_data.mass_spec.csv" - \
    --annotate "$(pwd)/curated_training_data.no_additional_ms.csv" - \
    --annotate "$(pwd)/t_cell_epitopes.csv" - \
    --fm-index-suffix .fm \
    --protein-column proteins_mouse

# Annotate viral proteins
time python annotate_proteins.py \
    "$(mhc2flurry-downloads path data_proteomes)/viral.uniprot.fasta.gz" \
    --annotate "$(pwd)/ms.by_pmid.csv" - \
    --annotate "$(pwd)/curated_training_data.csv" - \
    --annotate "$(pwd)/curated_training_data.mass_spec.csv" - \
    --annotate "$(pwd)/curated_training_data.no_additional_ms.csv" - \
    --annotate "$(pwd)/t_cell_epitopes.csv" - \
    --fm-index-suffix .fm \
    --protein-column proteins_viral


for i in $(ls *.csv)
do
    bzip2 $i
done

cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"
