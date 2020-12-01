import logging
logging.getLogger('tensorflow').disabled = True
logging.getLogger('matplotlib').disabled = True

import numpy
from numpy import testing
numpy.random.seed(0)
from tensorflow.random import set_seed
set_seed(2)

from nose.tools import eq_, assert_less, assert_greater, assert_almost_equal

import pandas
from sklearn.metrics import roc_auc_score

import mhcgnomes

from mhc2flurry.allele_encoding_pair import AlleleEncodingPair
from mhc2flurry.allele_encoding import AlleleEncoding
from mhc2flurry.class2_neural_network import Class2NeuralNetwork
from mhc2flurry.downloads import get_path
from mhc2flurry.common import random_peptides
from mhc2flurry.regression_target import from_ic50

from mhc2flurry.testing_utils import cleanup, startup
teardown = cleanup
setup = startup


# Fake pseudosequences
ALPHA_SEQUENCES = {
    "HLA-DRA*01:01": "AAAN",
}
BETA_SEQUENCES = {
    "HLA-DRB1*01:01": "AAAQ",
    "HLA-DRB1*03:01": "AAAK",
}


def make_allele_encoding_pair(allele_names):
    """
    Given a list of allele names, return an AlleleEncodingPair
    """
    parsed_alleles = pandas.Series([
        mhcgnomes.parse(name, infer_class2_pairing=True)
        for name in allele_names
    ])
    alpha = parsed_alleles.map(lambda p: p.alpha.to_string())
    beta = parsed_alleles.map(lambda p: p.beta.to_string())
    encoding = AlleleEncodingPair(
        AlleleEncoding(alpha, allele_to_sequence=ALPHA_SEQUENCES),
        AlleleEncoding(beta, allele_to_sequence=BETA_SEQUENCES),
    )
    return encoding


def test_memorize():
    train_df = pandas.DataFrame(
        {"peptide": random_peptides(200000, length=15)}
    ).set_index("peptide")
    train_df["HLA-DRB*01:01"] = 0
    train_df["HLA-DRB*03:01"] = 0

    # Each allele binds peptides matching some regex
    train_df.loc[
        train_df.index.str.contains("A.K"), "HLA-DRB*01:01"
    ] = 1.0
    train_df.loc[
        train_df.index.str.contains("Q.Q"), "HLA-DRB*03:01"
    ] = 1.0

    # Resample to have 1:1 binder / non-binder
    positive_train_df = train_df.loc[train_df.max(1) > 0.8]
    train_df = pandas.concat([
        positive_train_df,
        train_df.loc[~train_df.index.isin(positive_train_df.index)].sample(
            n=len(positive_train_df))
        ])

    print("Binders")
    print((train_df > 0.8).sum())

    print("Binder rate")
    print((train_df > 0.8).mean())

    stacked = train_df.stack().reset_index()
    stacked.columns = ['peptide', 'allele', 'measurement_value']

    # Memorize the dataset.
    allele_encoding = make_allele_encoding_pair(stacked.allele)
    model = Class2NeuralNetwork(
        random_negative_rate=0.0,  # setting this >0 seems to fail - why?
        layer_sizes=[8],
        allele_positionwise_embedding_size=2,
        patience=5,
        peptide_convolutions=[
            {'kernel_size': 3, 'filters': 8, 'activation': "relu"},
        ],
    )
    print(model.hyperparameters)

    model.fit(
        stacked.peptide.values,
        affinities=stacked["measurement_value"].values,
        allele_encoding_pair=allele_encoding
    )

    stacked["prediction"] = model.predict(
        stacked.peptide, allele_encoding_pair=allele_encoding)

    # Overall AUC
    stacked["binder"] = stacked.measurement_value > 0.8
    auc = roc_auc_score(stacked.binder, stacked.prediction)
    print("Overall AUC", auc)
    yield assert_greater, auc, 0.8

    # Can we discern a binder for one allele from another?
    binder_peptides = stacked.loc[stacked.binder].peptide.unique()
    stacked_binders = stacked.loc[stacked.peptide.isin(binder_peptides)]
    allele_specific_aucs = []
    for (allele, sub_df) in stacked_binders.groupby("allele"):
        print(allele)
        print(sub_df)
        auc = roc_auc_score(sub_df.binder, sub_df.prediction)
        allele_specific_aucs.append((allele, auc))
        yield assert_greater, auc, 0.8

    allele_specific_aucs = pandas.DataFrame(
        allele_specific_aucs, columns=["allele", "auc"])
    print("Allele specific AUCs:")
    print(allele_specific_aucs)

    #import ipdb ; ipdb.set_trace()

