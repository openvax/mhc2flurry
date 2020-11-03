"""
Filter and combine various peptide/MHC datasets to derive a composite training set,
optionally including eluted peptides identified by mass-spec.

The handle_pmid_XXXX functions should return a DataFrame with columns:
 - peptide
 - sample_id
 - hla [space separated list of alleles]
 - pulldown_antibody
 - format [monoallelic, multiallelic, DR-specific]
 - mhc_class [should be II]
 - sample type [an expression group, e.g. "spleen" or "expi293"]
 - cell_line [for samples deriving from a single known cell line]
"""
import sys
import argparse
import os
import json
import collections
from six.moves import StringIO

from mhcflurryii.common import normalize_allele_name


import pandas

parser = argparse.ArgumentParser(usage=__doc__)

parser.add_argument(
    "--ms-item",
    nargs="+",
    action="append",
    metavar="PMID FILE, ... FILE",
    default=[],
    help="Mass spec item to curate: PMID and list of files")
parser.add_argument(
    "--expression-item",
    nargs="+",
    action="append",
    metavar="LABEL FILE, ... FILE",
    default=[],
    help="Expression data to curate: dataset label and list of files")
parser.add_argument(
    "--ms-out",
    metavar="OUT.csv",
    help="Out file path (MS data)")
parser.add_argument(
    "--expression-out",
    metavar="OUT.csv",
    help="Out file path (RNA-seq expression)")
parser.add_argument(
    "--expression-metadata-out",
    metavar="OUT.csv",
    help="Out file path for expression metadata, i.e. which samples used")
parser.add_argument(
    "--debug",
    action="store_true",
    default=False,
    help="Leave user in pdb if PMID is unsupported")

PMID_HANDLERS = {}
EXPRESSION_HANDLERS = {}


def load(filenames, **kwargs):
    result = {}
    for filename in filenames:
        if filename.endswith(".csv"):
            result[filename] = pandas.read_csv(filename, **kwargs)
        elif filename.endswith(".xlsx") or filename.endswith(".xls"):
            result[filename] = pandas.read_excel(filename, **kwargs)
        else:
            result[filename] = filename

    return result


def debug(*filenames):
    loaded = load(filenames)
    import ipdb
    ipdb.set_trace()



PMID_31495665_SAMPLE_TYPES = {
        "HLA-DR_Lung": "lung",
        "HLA-DR_PBMC_HDSC": "pbmc",
        "HLA-DR_PBMC_RG1095": "pbmc",
        "HLA-DR_PBMC_RG1104": "pbmc",
        "HLA-DR_PBMC_RG1248": "pbmc",
        "HLA-DR_Spleen": "spleen",
        "MAPTAC_A*02:01": "mix:a375,expi293,hek293,hela",
        "MAPTAC_A*11:01": "mix:expi293,hela",
        "MAPTAC_A*32:01": "mix:a375,expi293,hela",
        "MAPTAC_B*07:02": "mix:a375,expi293,hela",
        "MAPTAC_B*45:01": "expi293",
        "MAPTAC_B*52:01": "mix:a375,expi293",
        "MAPTAC_C*03:03": "expi293",
        "MAPTAC_C*06:02": "mix:a375,expi293",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm+": "expi293",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm-": "expi293",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm+": "expi293",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm-": "expi293",
        "MAPTAC_DRB1*01:01": "mix:a375,b721,expi293,kg1,k562",
        "MAPTAC_DRB1*03:01": "expi293",
        "MAPTAC_DRB1*04:01": "expi293",
        "MAPTAC_DRB1*07:01": "mix:expi293,hek293",
        "MAPTAC_DRB1*11:01": "mix:expi293,k562,kg1",
        "MAPTAC_DRB1*12:01_dm+": "expi293",
        "MAPTAC_DRB1*12:01_dm-": "expi293",
        "MAPTAC_DRB1*15:01": "expi293",
        "MAPTAC_DRB3*01:01_dm+": "expi293",
        "MAPTAC_DRB3*01:01_dm-": "expi293",
}
CELL_LINE_MIXTURES = sorted(
    set(
        x for x in PMID_31495665_SAMPLE_TYPES.values()
        if x.startswith("mix:")))


def handle_pmid_25502872(filename):
    """Bergseng, ..., Sollid. Immunogenetics 2015 [PMID 25502872]"""
    return None


def handle_pmid_26495903(*filenames):
    """Sofron, ..., Fugmann. Eur. J. Immunol. 2015 [PMID 26495903]"""
    return None


def handle_pmid_26740625(*filenames):
    """Clement, ..., Santambrogio. J. Biol. Chem. 2016 [PMID 26740625]"""
    return None


def handle_pmid_27452731(*filenames):
    """Heyder, ..., Ytterberg. Mol. Cell. Proteomics 2016 [PMID 27452731]"""
    return None


def handle_pmid_27726376(*filenames):
    """Wang, ..., Costello. J. Proteom. Res. 2017"""
    return None


def handle_pmid_28329770(*filenames):
    """Khodadoust, ..., Alizadeh. Nature 2017 [PMID 28329770]"""
    return None


def handle_pmid_28467828(filename):
    """Ooi, ..., Kitching. Nature 2017 [PMID 28467828]"""
    return None


def handle_pmid_29314611(filename):
    """Ritz, ..., Fugmann. Proteomics 2018 [PMID 29314611]"""
    return None


def handle_pmid_29317506(*filenames):
    """Ting, ..., Rossjohn. J. Biol. Chem. 2018 [PMID 29317506]"""
    return None


def handle_pmid_29632711(*filenames):
    """Nelde, ..., Walz. Oncoimmunology 2018 [PMID 29632711]"""
    return None


def handle_pmid_31495665(filename):
    """Abelin, ..., Rooney Immunity 2019 [PMID 31495665]"""
    hla_type = {
        "HLA-DR_A375": None,
        "HLA-DR_Lung": "DRB1*01:01 DRB1*03:01 DRB3*01:01",
        "HLA-DR_PBMC_HDSC": "DRB1*03:01 DRB1*11:01 DRB3*01:01 DRB3*02:02",
        "HLA-DR_PBMC_RG1095": "HLA-DRA1*01:01-DRB1*03:01 HLA-DRA1*01:01-DRB1*11:01 HLA-DRA1*01:01-DRB3*01:01 HLA-DRA1*01:01-DRB3*02:02",
        "HLA-DR_PBMC_RG1104": "DRB1*01:01 DRB1*11:01 DRB3*02:02",
        "HLA-DR_PBMC_RG1248": "DRB1*03:01 DRB1*03:01 DRB3*01:01 DRB3*01:01",
        "HLA-DR_SILAC_Donor1_10minLysate": None,
        "HLA-DR_SILAC_Donor1_5hrLysate": None,
        "HLA-DR_SILAC_Donor1_DConly": None,
        "HLA-DR_SILAC_Donor1_UVovernight": None,
        "HLA-DR_SILAC_Donor2_DC_UV_16hr": None,
        "HLA-DR_SILAC_Donor2_DC_UV_24hr": None,
        "HLA-DR_Spleen": "DRB1*04:01 DRB4*01:03 DRB1*15:03 DRB5*01:01",
        "MAPTAC_A*02:01": "HLA-A*02:01",
        "MAPTAC_A*11:01": "HLA-A*11:01",
        "MAPTAC_A*32:01": "HLA-A*32:01",
        "MAPTAC_B*07:02": "HLA-B*07:02",
        "MAPTAC_B*45:01": "HLA-B*45:01",
        "MAPTAC_B*52:01": "HLA-B*52:01",
        "MAPTAC_C*03:03": "HLA-C*03:03",
        "MAPTAC_C*06:02": "HLA-C*06:02",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm+": "HLA-DPA1*01:03-DPB1*06:01",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm-": "HLA-DPA1*01:03-DPB1*06:01",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm+": "HLA-DQA1*01:02-DQB1*06:04",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm-": "HLA-DQA1*01:02-DQB1*06:04",
        "MAPTAC_DRB1*01:01": "HLA-DRA1*01:01-DRB1*01:01",
        "MAPTAC_DRB1*03:01": "HLA-DRA1*01:01-DRB1*03:01",
        "MAPTAC_DRB1*04:01": "HLA-DRA1*01:01-DRB1*04:01",
        "MAPTAC_DRB1*07:01": "HLA-DRA1*01:01-DRB1*07:01",
        "MAPTAC_DRB1*11:01": "HLA-DRA1*01:01-DRB1*11:01",
        "MAPTAC_DRB1*12:01_dm+": "HLA-DRA1*01:01-DRB1*12:01",
        "MAPTAC_DRB1*12:01_dm-": "HLA-DRA1*01:01-DRB1*12:01",
        "MAPTAC_DRB1*15:01": "HLA-DRA1*01:01-DRB1*15:01",
        "MAPTAC_DRB3*01:01_dm+": "HLA-DRA1*01:01-DRB3*01:01",
        "MAPTAC_DRB3*01:01_dm-": "HLA-DRA1*01:01-DRB3*01:01",
    }
    pulldown_antibody = {
        "HLA-DR_Lung": "L243 (HLA-DR)",
        "HLA-DR_PBMC_HDSC": "tal1b5 (HLA-DR)",
        "HLA-DR_PBMC_RG1095": "tal1b5 (HLA-DR)",
        "HLA-DR_PBMC_RG1104": "tal1b5 (HLA-DR)",
        "HLA-DR_PBMC_RG1248": "tal1b5 (HLA-DR)",
        "HLA-DR_Spleen": "L243 (HLA-DR)",
        "MAPTAC_A*02:01": "MAPTAC",
        "MAPTAC_A*11:01": "MAPTAC",
        "MAPTAC_A*32:01": "MAPTAC",
        "MAPTAC_B*07:02": "MAPTAC",
        "MAPTAC_B*45:01": "MAPTAC",
        "MAPTAC_B*52:01": "MAPTAC",
        "MAPTAC_C*03:03": "MAPTAC",
        "MAPTAC_C*06:02": "MAPTAC",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm+": "MAPTAC",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm-": "MAPTAC",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm+": "MAPTAC",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm-": "MAPTAC",
        "MAPTAC_DRB1*01:01": "MAPTAC",
        "MAPTAC_DRB1*03:01": "MAPTAC",
        "MAPTAC_DRB1*04:01": "MAPTAC",
        "MAPTAC_DRB1*07:01": "MAPTAC",
        "MAPTAC_DRB1*11:01": "MAPTAC",
        "MAPTAC_DRB1*12:01_dm+": "MAPTAC",
        "MAPTAC_DRB1*12:01_dm-": "MAPTAC",
        "MAPTAC_DRB1*15:01": "MAPTAC",
        "MAPTAC_DRB3*01:01_dm+": "MAPTAC",
        "MAPTAC_DRB3*01:01_dm-": "MAPTAC",
    }
    format = {
        "HLA-DR_Lung": "DR-specific",
        "HLA-DR_PBMC_HDSC": "DR-specific",
        "HLA-DR_PBMC_RG1095": "DR-specific",
        "HLA-DR_PBMC_RG1104": "DR-specific",
        "HLA-DR_PBMC_RG1248": "DR-specific",
        "HLA-DR_Spleen": "DR-specific",
        "MAPTAC_A*02:01": "monoallelic",
        "MAPTAC_A*11:01": "monoallelic",
        "MAPTAC_A*32:01": "monoallelic",
        "MAPTAC_B*07:02": "monoallelic",
        "MAPTAC_B*45:01": "monoallelic",
        "MAPTAC_B*52:01": "monoallelic",
        "MAPTAC_C*03:03": "monoallelic",
        "MAPTAC_C*06:02": "monoallelic",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm+": "monoallelic",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm-": "monoallelic",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm+": "monoallelic",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm-": "monoallelic",
        "MAPTAC_DRB1*01:01": "monoallelic",
        "MAPTAC_DRB1*03:01": "monoallelic",
        "MAPTAC_DRB1*04:01": "monoallelic",
        "MAPTAC_DRB1*07:01": "monoallelic",
        "MAPTAC_DRB1*11:01": "monoallelic",
        "MAPTAC_DRB1*12:01_dm+": "monoallelic",
        "MAPTAC_DRB1*12:01_dm-": "monoallelic",
        "MAPTAC_DRB1*15:01": "monoallelic",
        "MAPTAC_DRB3*01:01_dm+": "monoallelic",
        "MAPTAC_DRB3*01:01_dm-": "monoallelic",
    }
    mhc_class = {
        "HLA-DR_Lung": "II",
        "HLA-DR_PBMC_HDSC": "II",
        "HLA-DR_PBMC_RG1095": "II",
        "HLA-DR_PBMC_RG1104": "II",
        "HLA-DR_PBMC_RG1248": "II",
        "HLA-DR_Spleen": "II",
        "MAPTAC_A*02:01": "I",
        "MAPTAC_A*11:01": "I",
        "MAPTAC_A*32:01": "I",
        "MAPTAC_B*07:02": "I",
        "MAPTAC_B*45:01": "I",
        "MAPTAC_B*52:01": "I",
        "MAPTAC_C*03:03": "I",
        "MAPTAC_C*06:02": "I",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm+": "II",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm-": "II",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm+": "II",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm-": "II",
        "MAPTAC_DRB1*01:01": "II",
        "MAPTAC_DRB1*03:01": "II",
        "MAPTAC_DRB1*04:01": "II",
        "MAPTAC_DRB1*07:01": "II",
        "MAPTAC_DRB1*11:01": "II",
        "MAPTAC_DRB1*12:01_dm+": "II",
        "MAPTAC_DRB1*12:01_dm-": "II",
        "MAPTAC_DRB1*15:01": "II",
        "MAPTAC_DRB3*01:01_dm+": "II",
        "MAPTAC_DRB3*01:01_dm-": "II",
    }
    cell_line = {
        "HLA-DR_Lung": "",
        "HLA-DR_PBMC_HDSC": "",
        "HLA-DR_PBMC_RG1095": "",
        "HLA-DR_PBMC_RG1104": "",
        "HLA-DR_PBMC_RG1248": "",
        "HLA-DR_Spleen": "",
        "MAPTAC_A*02:01": "",
        "MAPTAC_A*11:01": "",
        "MAPTAC_A*32:01": "",
        "MAPTAC_B*07:02": "",
        "MAPTAC_B*45:01": "expi293",
        "MAPTAC_B*52:01": "",
        "MAPTAC_C*03:03": "expi293",
        "MAPTAC_C*06:02": "",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm+": "expi293",
        "MAPTAC_DPB1*06:01/DPA1*01:03_dm-": "expi293",
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm+": "expi293",  # don't actually see this in DataS1A!
        "MAPTAC_DQB1*06:04/DQA1*01:02_dm-": "expi293",
        "MAPTAC_DRB1*01:01": "",
        "MAPTAC_DRB1*03:01": "expi293",
        "MAPTAC_DRB1*04:01": "expi293",
        "MAPTAC_DRB1*07:01": "",
        "MAPTAC_DRB1*11:01": "",
        "MAPTAC_DRB1*12:01_dm+": "expi293",
        "MAPTAC_DRB1*12:01_dm-": "expi293",
        "MAPTAC_DRB1*15:01": "expi293",
        "MAPTAC_DRB3*01:01_dm+": "expi293",
        "MAPTAC_DRB3*01:01_dm-": "expi293",
    }

    df = pandas.read_excel(filename, sheet_name="DataS1B")
    results = []
    for sample_id in df.columns:
        if hla_type[sample_id] is None:
            print("Intentionally skipping", sample_id)
            continue

        result_df = pandas.DataFrame({
            "peptide": df[sample_id].dropna().values,
        })
        result_df["sample_id"] = sample_id
        result_df["hla"] = hla_type[sample_id]
        result_df["pulldown_antibody"] = pulldown_antibody[sample_id]
        result_df["format"] = format[sample_id]
        result_df["mhc_class"] = mhc_class[sample_id]
        result_df["sample_type"] = PMID_31495665_SAMPLE_TYPES[sample_id]
        result_df["cell_line"] = cell_line[sample_id]
        results.append(result_df)
    result_df = pandas.concat(results, ignore_index=True)
    result_df = result_df.loc[
        result_df.mhc_class == "II"
    ]
    return result_df


def handle_pmid_31611696(data_s1_filename, data_s2_filename):
    """Racle, ..., Gfeller. Nature Biotechnology 2019 [PMID 31611696]"""
    data_s1 = pandas.read_csv(data_s1_filename, sep=None).set_index("Sequence")
    data_s2 = pandas.read_csv(data_s2_filename, sep=None).set_index("Sequence")

    # HLA typing is given as a PDF in Supplementary Table 1.
    # In cases of ambiguous assignment we use the primary assignment.
    text = """
    3808_HMC MENINGIOMA DRB1*03:01 DRB1*07:01 DRB3*01:01 DRB4*01:01 DPA1*01:03 DPA1*02:01 DPB1*03:01 DPB1*11:01 DQA1*02:01 DQA1*05:01 DQB1*02:01 DQB1*02:02
    3830_NJF MENINGIOMA DRB1*04:04 DRB1*11:01 DRB3*02:02 DRB4*01:03 DPA1*01:03 DPB1*02:01 DPB1*06:01 DQA1*03:01 DQA1*05:05 DQB1*03:01 DQB1*03:02
    3849BR MENINGIOMA DRB1*11:04 DRB3*02:02 DPA1*01:03 DPB1*02:01 DPB1*04:01 DQA1*05:05 DQB1*03:01
    3865_DM MENINGIOMA DRB1*01:01 DRB1*07:01 DRB4*01:03 DPA1*01:03 DPB1*04:01 DPB1*20:01 DQA1*01:01 DQA1*02:01 DQB1*03:03 DQB1*05:01
    3869_GA MENINGIOMA DRB1*01:03 DRB1*04:04 DRB4*01:03 DPA1*01:03 DPB1*04:01 DPB1*126:01 DQA1*03:01 DQA1*05:05 DQB1*03:01 DQB1*03:02
    3911_ME MENINGIOMA DRB1*11:01 DRB3*02:02 DPA1*01:03 DPB1*04:01 DQA1*05:05 DQB1*03:01
    3912_BAM MENINGIOMA DRB1*03:01 DRB1*04:01 DRB3*01:01 DRB4*01:03 DPA1*01:03 DPB1*04:01 DQA1*03:01 DQA1*05:01 DQB1*02:01 DQB1*03:02
    3947_GA MENINGIOMA DRB1*01:01 DRB1*13:01 DRB3*01:01 DPA1*01:03 DPB1*02:01 DPB1*04:02 DQA1*01:01 DQA1*01:03 DQB1*05:01 DQB1*06:03
    3971_ORA MENINGIOMA DRB1*13:03 DRB1*07:01 DRB3*01:01 DRB4*01:01 DPA1*01:03 DPA1*02:02 DPB1*04:01 DQA1*02:01 DQA1*05:05 DQB1*02:02 DQB1*03:01
    3993 MENINGIOMA DRB1*07:01 DRB1*15:01 DRB4*01:03 DRB5*01:01 DPA1*01:03 DPA1*02:01 DPB1*04:01 DPB1*17:01 DQA1*01:02 DQA1*02:01 DQB1*02:02 DQB1*06:02
    4001 MENINGIOMA DRB1*13:01 DRB1*14:01 DRB3*01:01 DRB3*02:02 DPA1*01:03 DPB1*04:01 DPB1*04:02 DQA1*01:03 DQA1*01:04 DQB1*05:03 DQB1*06:03
    4021 MENINGIOMA DRB1*11:01 DRB1*04:05 DRB3*02:02 DRB4*01:03 DPA1*01:03 DPB1*03:01 DPB1*104:01 DQA1*03:03 DQA1*05:05 DQB1*02:02 DQB1*03:01
    4037_DC MENINGIOMA DRB1*01:01 DPA1*01:03 DPB1*04:01 DPB1*06:01 DQA1*01:01 DQB1*05:01
    4052_BA MENINGIOMA DRB1*03:01 DRB1*11:04 DRB3*01:01 DRB3*02:02 DPA1*01:03 DPB1*04:01 DQA1*05:01 DQA1*05:05 DQB1*02:01 DQB1*03:01
    BP455 B-CELL DRB1*10:01 DRB1*13:01 DRB3*01:01 DPA1*01:03 DPB1*02:01 DQA1*01:05 DQA1*01:10 DQB1*05:01 DQB1*06:03
    CD165 B-CELL DRB1*11:01 DRB3*02:02 DPA1*01:03 DPB1*04:01 DPB1*04:02 DQA1*05:05 DQB1*03:01
    CM647 B-CELL DRB1*07:01 DRB1*16:01 DRB4*01:03 DRB5*02:02 DPA1*01:03 DPB1*02:01 DPB1*23:01 DQA1*01:02 DQA1*02:01 DQB1*02:02 DQB1*05:02
    GD149 B-CELL DRB1*07:01 DRB1*13:01 DRB3*01:01 DRB4*01:01 DPA1*01:03 DPA1*02:01 DPB1*03:01 DPB1*04:01 DQA1*01:10 DQA1*02:01 DQB1*02:02 DQB1*06:03
    JY B-CELL DRB1*04:04 DRB1*13:01 DRB3*01:01 DRB4*01:03 DPA1*01:03 DPB1*02:01 DPB1*04:01 DQA1*01:03 DQA1*03:01 DQB1*03:02 DQB1*06:03
    PD42 B-CELL DRB1*01:02 DRB1*15:01 DRB5*01:01 DPA1*01:03 DPA1*02:02 DPB1*04:01 DPB1*05:01 DQA1*01:01 DQA1*01:02 DQB1*05:01 DQB1*06:02
    RA957 B-CELL DRB1*04:01 DRB1*08:01 DRB4*01:03 DPA1*01:03 DPB1*04:01 DPB1*04:02 DQA1*03:03 DQA1*04:01 DQB1*03:01 DQB1*04:02
    TIL1 TIL DRB1*01:01 DRB1*04:08 DRB4*01:03 DPA1*01:03 DPB1*02:01 DPB1*04:01 DQA1*01:01 DQA1*03:03 DQB1*03:01 DQB1*05:01
    TIL3 TIL DRB1*12:01 DRB1*15:01 DRB3*02:02 DRB5*01:01 DPA1*01:03 DPB1*03:01 DPB1*04:01 DQA1*01:02 DQA1*05:05 DQB1*03:01 DQB1*05:02
    """
    rows = [
        row.split() for row in text.strip().split("\n")
    ]
    rows = [
        (row[0].replace("_", "-"), row[1], " ".join(row[2:])) for row in rows
    ]
    info_df = pandas.DataFrame(rows, columns=["kind", "sample_type", "hla"])
    info_df = info_df.set_index("kind")

    # Data S1
    renames = {
        c : c.replace("Intensity", "").replace("_II", "").strip()
        for c in data_s1.columns if c.startswith("Intensity")
    }

    data_s1 = data_s1[sorted(renames)].rename(columns=renames).rename(columns={
        "3830NJF": "3830-NJF",
        "3865DM": "3865-DM",
        "3912BAM": "3912-BAM",
        "3865DM": "3865-DM",
        "CD165_ IFNg": "CD165_IFNg",
    })

    result1_df = data_s1.stack().reset_index()
    result1_df.columns = ["peptide", "sample_id", "intensity"]
    result1_df = result1_df.loc[result1_df.intensity > 0]
    result1_df["kind"] = result1_df.sample_id.map(lambda s: {
        "JY_DR": "JY",
        "CD165_IFNg": "CD165",
    }.get(s, s))
    result1_df["hla"] = result1_df.kind.map(info_df.hla)
    result1_df["pulldown_antibody"] = "HB145"
    result1_df["format"] = "MULTIALLELIC"
    result1_df.loc[
        result1_df.sample_id == "JY_DR",
        "format"
    ] = "DR-specific"
    result1_df["mhc_class"] = "II"
    result1_df["sample_type"] = result1_df.kind.map(info_df.sample_type)
    result1_df["cell_line"] = [
        row.kind if row.sample_type == "B-CELL" else ""
        for _, row in result1_df.iterrows()
    ]
    del result1_df["kind"]

    # Data S2
    renames = {
        c : c.replace("Intensity", "").replace("_II", "").strip()
        for c in data_s2.columns if c.startswith("Intensity")
    }

    data_s2 = data_s2[sorted(renames)].rename(columns=renames).rename(columns={
        "3830NJF": "3830-NJF",
        "3865DM": "3865-DM",
        "3912BAM": "3912-BAM",
        "3865DM": "3865-DM",
        "CD165_ IFNg": "CD165_IFNg",
    })
    result2_df = data_s2.stack().reset_index()
    result2_df.columns = ["peptide", "sample_id", "intensity"]
    result2_df["kind"] = result2_df.sample_id.str.replace(
        "-HLA-DR", "").str.replace("-depleted", "").str.replace("_", "-")
    result2_df["hla"] = result2_df.kind.map(info_df.hla)
    result2_df["pulldown_antibody"] = ""
    assert all(result2_df.sample_id.map(
        lambda s: s.endswith("DR-depleted") or s.endswith("-DR")))
    result2_df["format"] = result2_df.sample_id.map(
        lambda s: "DR-depleted" if "DR-depleted" in s else "DR-specific")

    result2_df["mhc_class"] = "II"
    result2_df["sample_type"] = result2_df.kind.map(info_df.sample_type)
    result2_df["cell_line"] = [
        row.kind if row.sample_type == "B-CELL" else "" for _, row in
        result2_df.iterrows()
    ]
    del result2_df["kind"]
    result_df = pandas.concat([result1_df, result2_df], ignore_index=True)

    # DR-specific samples used HB298 antibody
    result_df.loc[
        result_df.format == "DR-specific",
        "pulldown_antibody"
    ] = "HB298"

    # Subsample alleles to just DR alleles for DR-specific samples.
    result_df.loc[
        result_df.format == "DR-specific",
        "hla"
    ] = result_df.loc[result_df.format == "DR-specific", "hla"].map(
        lambda s: " ".join([allele for allele in s.split() if "DR" in allele])
    )

    return result_df


def Xhandle_pmid_27869121(filename):
    """Bassani-Sternberg, ..., Krackhardt Nature Comm. 2016 [PMID 27869121]"""
    # Although this dataset has class II data also, we are only extracting
    # class I for now.
    df = pandas.read_excel(filename, skiprows=1)

    # Taking these from:
    # Supplementary Table 2: Information of patients selected for neoepitope
    # identification
    # For the Mel5 sample, only two-digit alleles are shown (A*01, A*25,
    # B*08, B*18) so we are skipping that sample for now.
    hla_df = pandas.DataFrame([
        ("Mel-8", "HLA-A*01:01 HLA-A*03:01 HLA-B*07:02 HLA-B*08:01 HLA-C*07:01 HLA-C*07:02"),
        ("Mel-12", "HLA-A*01:01 HLA-B*08:01 HLA-C*07:01"),
        ("Mel-15", "HLA-A*03:01 HLA-A*68:01 HLA-B*27:05 HLA-B*35:03 HLA-C*02:02 HLA-C*04:01"),
        ("Mel-16", "HLA-A*01:01 HLA-A*24:02 HLA-B*07:02 HLA-B*08:01 HLA-C*07:01 HLA-C*07:02"),
    ], columns=["sample_id", "hla"]).set_index("sample_id")

    # We assert below that none of the class I hit peptides were found in any
    # of the class II pull downs.
    class_ii_cols = [
        c for c in df.columns if c.endswith("HLA-II (arbitrary units)")
    ]
    class_ii_hits = set(df.loc[
        (df[class_ii_cols].fillna(0.0).sum(1) > 0)
    ].Sequence.unique())

    results = []
    for (sample_id, hla) in hla_df.hla.items():
        intensity_col = "Intensity %s_HLA-I (arbitrary units)" % sample_id
        sub_df = df.loc[
            (df[intensity_col].fillna(0.0) > 0)
        ]
        filtered_sub_df = sub_df.loc[
            (~sub_df.Sequence.isin(class_ii_hits))
        ]
        peptides = filtered_sub_df.Sequence.unique()
        assert not any(p in class_ii_hits for p in peptides)

        result_df = pandas.DataFrame({
            "peptide": peptides,
        })
        result_df["sample_id"] = sample_id
        result_df["hla"] = hla_df.loc[sample_id, "hla"]
        result_df["pulldown_antibody"] = "W6/32"
        result_df["format"] = "multiallelic"
        result_df["mhc_class"] = "I"
        result_df["sample_type"] = "melanoma_met"
        result_df["cell_line"] = None
        results.append(result_df)

    result_df = pandas.concat(results, ignore_index=True)
    return result_df




EXPRESSION_GROUPS_ROWS = []


def make_expression_groups(dataset_identifier, df, groups):
    result_df = pandas.DataFrame(index=df.index)
    for (label, columns) in groups.items():
        for col in columns:
            if col not in df.columns:
                raise ValueError(
                    "Missing: %s. Available: %s" % (col, df.columns.tolist()))
        result_df[label] = df[columns].mean(1)
        EXPRESSION_GROUPS_ROWS.append((dataset_identifier, label, columns))
    return result_df


def handle_expression_GSE113126(*filenames):
    """
    Barry, ..., Krummel Nature Medicine 2018 [PMID 29942093]

    This is the melanoma met RNA-seq dataset.

    """

    df = pandas.read_csv(filenames[0], sep="\t", index_col=0)
    df = df[[]]  # no columns

    for filename in filenames:
        df[os.path.basename(filename)] = pandas.read_csv(
            filename, sep="\t", index_col=0)["TPM"]

    assert len(df.columns) == len(filenames)

    groups = {
        "sample_type:MELANOMA_MET": df.columns.tolist(),
    }
    return [make_expression_groups("GSE113126", df, groups)]


def handle_expression_expression_atlas_22460905(filename):
    df = pandas.read_csv(filename, sep="\t", skiprows=4, index_col=0)
    del df["Gene Name"]
    df.columns = df.columns.str.lower()
    df = df.fillna(0.0)

    def matches(*strings):
        return [c for c in df.columns if all(s in c for s in strings)]

    groups = {
        "sample_type:B-LCL": (
            matches("b-cell", "lymphoblast") + matches("b acute lymphoblastic")),
        "sample_type:B-CELL": matches("b-cell"),
        "sample_type:B721-LIKE": matches("b-cell"),
        "sample_type:MELANOMA_CELL_LINE": matches("melanoma"),
        "sample_type:MELANOMA": matches("melanoma"),
        "sample_type:A375-LIKE": matches("melanoma"),
        "sample_type:KG1-LIKE": matches("myeloid leukemia"),

        # Using a fibrosarcoma cell line for our fibroblast sample.
        "sample_type:FIBROBLAST": ['fibrosarcoma, ht-1080'],

        # For GBM tissue we are just using a mixture of cell lines.
        "sample_type:GLIOBLASTOMA_TISSUE": matches("glioblastoma"),

        "cell_line:THP-1": ["childhood acute monocytic leukemia, thp-1"],
        "cell_line:HL-60": ["adult acute myeloid leukemia, hl-60"],
        "cell_line:U-87": ['glioblastoma, u-87 mg'],
        "cell_line:LNT-229": ['glioblastoma, ln-229'],
        "cell_line:T98G": ['glioblastoma, t98g'],
        "cell_line:SK-MEL-5": ['cutaneous melanoma, sk-mel-5'],
        'cell_line:MEWO': ['melanoma, mewo'],
        "cell_line:HCC1937": ['breast ductal adenocarcinoma, hcc1937'],
        "cell_line:HCT116": ['colon carcinoma, hct 116'],
        "cell_line:HCC1143": ['breast ductal adenocarcinoma, hcc1143'],
    }
    return [make_expression_groups("expression_atlas_22460905", df, groups)]


def handle_expression_human_protein_atlas(*filenames):
    (cell_line_filename,) = [f for f in filenames if "celline" in f]
    (blood_filename,) = [f for f in filenames if "blood" in f]
    (gtex_filename,) = [f for f in filenames if "gtex" in f]

    cell_line_df = pandas.read_csv(cell_line_filename, sep="\t")
    blood_df = pandas.read_csv(blood_filename, sep="\t", index_col=0)
    gtex_df = pandas.read_csv(gtex_filename, sep="\t")

    cell_line_df = cell_line_df.pivot(
        index="Gene", columns="Cell line", values="TPM")

    gtex_df = gtex_df.pivot(
        index="Gene", columns="Tissue", values="TPM")

    return [
        make_expression_groups(
            "human_protein_atlas:%s" % os.path.basename(blood_filename),
            blood_df,
            groups={
                "sample_type:PBMC": [
                    c for c in blood_df.columns if "total PBMC" in c
                ],

                # for samples labeled leukapheresis we also use PBMC
                "sample_type:LEUKAPHERESIS": [
                    c for c in blood_df.columns if "total PBMC" in c
                ],

                # for samples labeled TIL we are also using PBMC
                "sample_type:TIL": [
                    c for c in blood_df.columns if "total PBMC" in c
                ],
            }),
        make_expression_groups(
            "human_protein_atlas:%s" % os.path.basename(cell_line_filename),
            cell_line_df,
            groups={
                "cell_line:HELA": ['HeLa'],
                "cell_line:K562": ["K-562"],
                "cell_line:HEK293": ['HEK 293'],
                "cell_line:RPMI8226": ['RPMI-8226'],
                "cell_line:EXPI293": ['HEK 293'],  # EXPI293 derived from HEK293
            }),
        make_expression_groups(
            "human_protein_atlas:%s" % os.path.basename(gtex_filename),
            gtex_df,
            groups={
                "sample_type:LUNG": ["lung"],
                "sample_type:SPLEEN": ["spleen"],
                "sample_type:OVARY": ["ovary"],
                "sample_type:KIDNEY": ["kidney"],

                # This is bad! I just can't find anything better currently.
                # We should find some meningioma RNA-seq and switch to that.
                "sample_type:MENINGIOMA": [
                    "amygdala", "basal ganglia", "cerebellum", "cerebral cortex",
                    "midbrain", "spinal cord",
                ],
            }),
    ]


def make_expression_mixtures(expression_df):
    global CELL_LINE_MIXTURES
    groups = {}
    for mix in CELL_LINE_MIXTURES:
        components = []
        for item in mix.replace("mix:", "").upper().split(","):
            if "cell_line:%s" % item in expression_df.columns:
                components.append("cell_line:%s" % item)
            else:
                print("No cell line, falling back on similar: ", item)
                components.append("sample_type:%s-LIKE" % item)
        groups["sample_type:" + mix.upper()] = components
    missing = set()
    for some in groups.values():
        for item in some:
            if item not in expression_df.columns:
                missing.add(item)
    if missing:
        raise ValueError(
            "Missing [%d]: %s. Available: %s" % (
                len(missing), missing, expression_df.columns.tolist()))
    return make_expression_groups("mixtures", expression_df, groups)


# Add all functions with names like handle_pmid_XXXX to PMID_HANDLERS dict.
for (key, value) in list(locals().items()):
    if key.startswith("handle_pmid_"):
        PMID_HANDLERS[key.replace("handle_pmid_", "")] = value
    elif key.startswith("handle_expression_"):
        EXPRESSION_HANDLERS[key.replace("handle_expression_", "")] = value


def run():
    args = parser.parse_args(sys.argv[1:])

    expression_dfs = []
    for (i, item_tpl) in enumerate(args.expression_item):
        (label, filenames) = (item_tpl[0], item_tpl[1:])
        label = label.replace("-", "_")
        print(
            "Processing expression item %d of %d" % (i + 1, len(args.expression_item)),
            label,
            *[os.path.abspath(f) for f in filenames])

        expression_dfs_for_item = []
        handler = None
        if label in EXPRESSION_HANDLERS:
            handler = EXPRESSION_HANDLERS[label]
            expression_dfs_for_item = handler(*filenames)
        elif args.debug:
            debug(*filenames)
        else:
            raise NotImplementedError(label)

        if expression_dfs_for_item:
            print(
                "Processed expression data",
                label,
                "result dataframes",
                len(expression_dfs_for_item))
            print(*[e.columns for e in expression_dfs_for_item])
            expression_dfs.extend(expression_dfs_for_item)

    expression_df = expression_dfs[0]
    for other in expression_dfs[1:]:
        expression_df = pandas.merge(
            expression_df, other, how='outer', left_index=True, right_index=True)

    print("Genes in each expression dataframe: ",
        *[len(e) for e in expression_dfs])
    print("Genes in merged expression dataframe", len(expression_df))

    if CELL_LINE_MIXTURES:
        print("Generating cell line mixtures.")
        expression_mixture_df = make_expression_mixtures(expression_df)
        expression_df = pandas.merge(
            expression_df,
            expression_mixture_df,
            how='outer',
            left_index=True,
            right_index=True)

    ms_dfs = []
    for (i, item_tpl) in enumerate(args.ms_item):
        (pmid, filenames) = (item_tpl[0], item_tpl[1:])
        print(
            "Processing MS item %d of %d" % (i + 1, len(args.ms_item)),
            pmid,
            *[os.path.abspath(f) for f in filenames])

        ms_df = None
        handler = None
        if pmid in PMID_HANDLERS:
            handler = PMID_HANDLERS[pmid]
            ms_df = handler(*filenames)
        elif args.debug:
            debug(*filenames)
        else:
            raise NotImplementedError(pmid)

        if ms_df is not None:
            ms_df["pmid"] = pmid
            if "original_pmid" not in ms_df.columns:
                ms_df["original_pmid"] = pmid
            if "expression_dataset" not in ms_df.columns:
                ms_df["expression_dataset"] = ""
            ms_df = ms_df.applymap(str).applymap(str.upper)
            ms_df["sample_id"] = ms_df.sample_id.str.replace(" ", "")
            print("*** PMID %s: %d peptides ***" % (pmid, len(ms_df)))
            if handler is not None:
                print(handler.__doc__)
            print("Counts by sample id:")
            print(ms_df.groupby("sample_id").peptide.nunique())
            print("")
            print("Counts by sample type:")
            print(ms_df.groupby("sample_type").peptide.nunique())
            print("****************************")

            for value in ms_df.expression_dataset.unique():
                if value and value not in expression_df.columns:
                    raise ValueError("No such expression dataset", value)

            ms_dfs.append(ms_df)
        else:
            print("Skipping MS item", pmid)

    ms_df = pandas.concat(ms_dfs, ignore_index=True, sort=False)
    ms_df["cell_line"] = ms_df["cell_line"].fillna("")
    ms_df["hla"] = ms_df["hla"].str.strip().str.replace(r'\s+', ' ').map(
        lambda hla: " ".join(
            [
                normalize_allele_name(a, raise_on_error=True)
                for a in hla.split()
            ]))

    sample_table = ms_df[
        ["sample_id", "pmid", "expression_dataset", "cell_line", "sample_type"]
    ].drop_duplicates().set_index("sample_id")

    sample_id_to_expression_dataset = sample_table.expression_dataset.to_dict()
    for (sample_id, value) in sorted(sample_id_to_expression_dataset.items()):
        if value:
            print("Expression dataset for sample", sample_id, "already assigned")
            continue
        cell_line_col = "cell_line:" + sample_table.loc[sample_id, "cell_line"]
        sample_type_col = "sample_type:" + (
            sample_table.loc[sample_id, "sample_type"])

        expression_dataset = None
        for col in [cell_line_col, sample_type_col]:
            if col in expression_df.columns:
                expression_dataset = col
                break

        if not expression_dataset:
            print("*" * 20)
            print("No expression dataset for sample ", sample_id)
            print("Sample info:")
            print(sample_table.loc[sample_id])
            print("*" * 20)

        sample_id_to_expression_dataset[sample_id] = expression_dataset
        print(
            "Sample", sample_id, "assigned exp. dataset", expression_dataset)

    print("Expression dataset usage:")
    print(pandas.Series(sample_id_to_expression_dataset).value_counts())

    missing = [
        key for (key, value) in
        sample_id_to_expression_dataset.items()
        if value is None
    ]
    if missing:
        print("Missing expression data for samples", *missing)
        print(
            "Missing cell lines: ",
            *sample_table.loc[missing, "cell_line"].dropna().drop_duplicates().tolist())
        print("Missing sample types: ", *sample_table.loc[
            missing, "sample_type"].dropna().drop_duplicates().tolist())
        if args.debug:
            import ipdb; ipdb.set_trace()
        else:
            raise ValueError("Missing expression data for samples: ", missing)

    ms_df["expression_dataset"] = ms_df.sample_id.map(
        sample_id_to_expression_dataset)

    cols = [
        "pmid",
        "sample_id",
        "peptide",
        "format",
        "mhc_class",
        "hla",
        "expression_dataset",
    ]
    cols += [c for c in sorted(ms_df.columns) if c not in cols]
    ms_df = ms_df[cols]

    null_df = ms_df.loc[ms_df.isnull().any(1)]
    if len(null_df) > 0:
        print("Nulls:")
        print(null_df)
    else:
        print("No nulls.")

    # Each sample should be coming from only one experiment.
    assert ms_df.groupby("sample_id").pmid.nunique().max() == 1, (
        ms_df.groupby("sample_id").pmid.nunique().sort_values())

    expression_df.to_csv(args.expression_out, index=True)
    print("Wrote: %s" % os.path.abspath(args.expression_out))

    ms_df.to_csv(args.ms_out, index=False)
    print("Wrote: %s" % os.path.abspath(args.ms_out))

    if args.expression_metadata_out is not None:
        expression_metadata_df = pandas.DataFrame(
            EXPRESSION_GROUPS_ROWS,
            columns=["expression_dataset", "label", "samples"])
        expression_metadata_df["samples"] = expression_metadata_df[
            "samples"
        ].map(json.dumps)
        expression_metadata_df.to_csv(args.expression_metadata_out, index=False)
        print("Wrote: %s" % os.path.abspath(args.expression_metadata_out))


if __name__ == '__main__':
    run()
