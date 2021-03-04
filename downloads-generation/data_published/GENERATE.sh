#!/bin/bash
#
# Download published MHC II ligand data.
#
set -e
set -x

DOWNLOAD_NAME=data_published
SCRATCH_DIR=${TMPDIR-/tmp}/mhc2flurry-downloads-generation
SCRIPT_ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR=$(dirname "$SCRIPT_ABSOLUTE_PATH")

mkdir -p "$SCRATCH_DIR"
rm -rf "$SCRATCH_DIR/$DOWNLOAD_NAME"
mkdir "$SCRATCH_DIR/$DOWNLOAD_NAME"

# Send stdout and stderr to a logfile included with the archive.
exec >  >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt")
exec 2> >(tee -ia "$SCRATCH_DIR/$DOWNLOAD_NAME/LOG.txt" >&2)

date

cd $SCRATCH_DIR/$DOWNLOAD_NAME

############################################
# MS: Class II
############################################

## Scally, ..., Rossjohn. J. Exp. Med. [PMID 24190431]
# PMID=24190431
# mkdir -p ms/$PMID
# TODO: couldn't find peptides, need to look again.

# Bergseng, ..., Sollid. Immunogenetics 2015 [PMID 25502872]
PMID=25502872
mkdir -p ms/$PMID
wget -q 'https://static-content.springer.com/esm/art%3A10.1007%2Fs00251-014-0819-9/MediaObjects/251_2014_819_MOESM3_ESM.xlsx' -P ms/$PMID

# Sofron, ..., Fugmann. Eur. J. Immunol. 2015 [PMID 26495903]
PMID=26495903
mkdir -p ms/$PMID
wget -q 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Feji.201545930&file=eji3483-sup-0003-supinfo.xlsx' -P ms/$PMID
wget -q 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Feji.201545930&file=eji3483-sup-0004-supinfo.xlsx' -P ms/$PMID
wget -q 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Feji.201545930&file=eji3483-sup-0005-supinfo.xlsx' -P ms/$PMID

# Bassani-Sternberg, ..., Krackhardt Nature Comm. 2016 [PMID 27869121]
PMID=27869121
mkdir -p ms/$PMID
wget -q "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms13404/MediaObjects/41467_2016_BFncomms13404_MOESM1318_ESM.xlsx" -P ms/$PMID

# Clement, ..., Santambrogio. J. Biol. Chem. 2016 [PMID 26740625]
PMID=26740625
mkdir -p ms/$PMID
wget -q 'https://www.jbc.org/cms/10.1074/jbc.M115.655738/attachment/687abf48-576f-41e1-8f9f-dc40fbfa27ae/mmc1.zip' -P ms/$PMID
pushd ms/$PMID
unzip *.zip
rm *.zip
popd

# Heyder, ..., Ytterberg. Mol. Cell. Proteomics 2016 [PMID 27452731]
PMID=27452731
mkdir -p ms/$PMID
wget -q 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013314/bin/10.1074_M116.060764_mcp.M116.060764-5.xlsx' -P ms/$PMID
wget -q 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013314/bin/10.1074_M116.060764_mcp.M116.060764-6.xlsx' -P ms/$PMID
wget -q 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013314/bin/10.1074_M116.060764_mcp.M116.060764-7.xlsx' -P ms/$PMID
wget -q 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013314/bin/10.1074_M116.060764_mcp.M116.060764-8.xlsx' -P ms/$PMID
wget -q 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013314/bin/10.1074_M116.060764_mcp.M116.060764-9.xlsx' -P ms/$PMID
wget -q 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5013314/bin/10.1074_M116.060764_mcp.M116.060764-5.xlsx' -P ms/$PMID

## Wang, ..., Costello. J. Proteom. Res. 2017 [PMID 27726376]
PMID=27726376
mkdir -p ms/$PMID
wget -q 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00386/suppl_file/pr6b00386_si_001.docx' -P ms/$PMID
wget -q 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00386/suppl_file/pr6b00386_si_002.xlsx' -P ms/$PMID
wget -q 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00386/suppl_file/pr6b00386_si_003.xlsx' -P ms/$PMID
wget -q 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00386/suppl_file/pr6b00386_si_004.xlsx' -P ms/$PMID

## Khodadoust, ..., Alizadeh. Nature 2017 [PMID 28329770]
# PMID=28329770
# mkdir -p ms/$PMID
# TODO. PRIDE Archive at www.ebi.ac.uk/pride/archive under accession numbers PXD004746 and PXD005704

# Ooi, ..., Kitching. Nature 2017 [PMID 28467828]
PMID=28467828
mkdir -p ms/$PMID
wget -q 'https://static-content.springer.com/esm/art%3A10.1038%2Fnature22329/MediaObjects/41586_2017_BFnature22329_MOESM1_ESM.xlsx' -P ms/$PMID

# Ritz, ..., Fugmann. Proteomics 2018 [PMID 29314611]
PMID=29314611
mkdir -p ms/$PMID
wget -q 'https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpmic.201700246&file=pmic12799-sup-0006-TableS4.xlsx' -P ms/$PMID

# Ting, ..., Rossjohn. J. Biol. Chem. 2018 [PMID 29317506]
PMID=29317506
mkdir -p ms/$PMID
wget -q 'https://www.jbc.org/cms/10.1074/jbc.RA117.001013/attachment/3a33375f-7acd-420c-ac79-e3adbb6cf394/mmc1.zip' -P ms/$PMID
pushd ms/$PMID
unzip *.zip
rm *.zip
popd

# Nelde, ..., Walz. Oncoimmunology 2018 [PMID 29632711]
PMID=29632711
mkdir -p ms/$PMID
wget -q 'https://www.tandfonline.com/doi/suppl/10.1080/2162402X.2017.1316438/suppl_file/koni_a_1316438_sm6974.zip' -P ms/$PMID
pushd ms/$PMID
unzip *.zip
rm *.zip
popd

## Alvaro-Benito, ..., Freund. Front. Immunol 2018 [PMID 29774024]
# PMID=29774024
# mkdir -p ms/$PMID
# TODO

## Nanaware, ..., Stern. Mol. Cell. Proteomics 2019 [PMID 30573663]
# PMID=30573663
# mkdir -p ms/$PMID
# TODO: Accession MSV000082570 at http://massive.ucsd.edu

# Abelin, ..., Rooney. Immunity 2019 [PMID 31495665]
PMID=31495665
mkdir -p ms/$PMID
wget -q 'https://ars.els-cdn.com/content/image/1-s2.0-S1074761319303632-mmc2.xlsx' -P ms/$PMID

# Racle, ..., Gfeller. Nature Biotechnology 2019 [PMID 31611696]
PMID=31611696
mkdir -p ms/$PMID
wget -q 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0289-6/MediaObjects/41587_2019_289_MOESM4_ESM.txt' -P ms/$PMID
wget -q 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0289-6/MediaObjects/41587_2019_289_MOESM5_ESM.txt' -P ms/$PMID

############################################
# Non MS: T cell epitopes, etc.
############################################
#
# Reynissonm ..., Nielsen. Nucleic Acids Res. 2020 [PMID 32406916]
#
PMID=32406916
mkdir -p other/$PMID
wget -q 'https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/suppl/CD4_epitopes/CD4_epitopes.fsa' -P other/$PMID
wget -q 'http://www.cbs.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/NetMHCIIpan_train.tar.gz' -P other/$PMID
pushd other/$PMID
tar xvzf *.tar.gz
rm *.tar.gz
for i in $(find . -mindepth 1 -type f)
do
    bzip2 $i
done
popd

############################################
# RNA-seq expression data (TPMs)
############################################
# CCLE as processed by expression atlas
DATASET=expression-atlas-22460905
mkdir -p expression/$DATASET
wget -q https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2770/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv -P expression/$DATASET

# Human protein atlas
DATASET=human-protein-atlas
mkdir -p expression/$DATASET
pushd expression/$DATASET
wget -q https://www.proteinatlas.org/download/rna_celline.tsv.zip
wget -q https://www.proteinatlas.org/download/rna_blood_cell_sample_tpm_m.tsv.zip
wget -q https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip
for i in $(ls *.zip)
do
    unzip $i
    rm $i
done
popd

# Melanoma. Original publication
# Barry, ..., Krummel Nature Medicine 2018 [PMID 29942093].
DATASET=GSE113126
mkdir -p expression/$DATASET 
pushd expression/$DATASET
wget -q "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE113126&format=file" -O GSE113126_RAW.tar
tar -xvf GSE113126_RAW.tar
rm GSE113126_RAW.tar
popd

############################################
cp $SCRIPT_ABSOLUTE_PATH .
bzip2 LOG.txt
RESULT="$SCRATCH_DIR/${DOWNLOAD_NAME}.$(date +%Y%m%d).tar.bz2"
tar -cjf "$RESULT" *
echo "Created archive: $RESULT"

