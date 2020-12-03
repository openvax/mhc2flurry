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


def make_allele_encoding_pair(allele_names, alpha_sequences, beta_sequences):
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
        AlleleEncoding(alpha, allele_to_sequence=alpha_sequences),
        AlleleEncoding(beta, allele_to_sequence=beta_sequences),
    )
    return encoding


def Xtest_simple():
    # Fake pseudosequences
    alpha_sequences = {
        "HLA-DRA*01:01": "AAAN",
    }
    beta_sequences = {
        "HLA-DRB1*01:01": "AAAQ",
        "HLA-DRB1*03:01": "AAAK",
    }
    motifs = {
        "HLA-DRB1*01:01": "A.K",
        "HLA-DRB1*03:01": "Q.Q",
    }

    df = pandas.DataFrame(
        {"peptide": random_peptides(200000, length=15)}
    ).set_index("peptide")

    for (allele, motif) in motifs.items():
        df[allele] = (df.index.str.contains(motif)).astype(int)

    # Resample to have 1:1 binder / non-binder
    positive_train_df = df.loc[df.max(1) > 0.8]
    train_df = pandas.concat([
        positive_train_df,
        df.loc[~df.index.isin(positive_train_df.index)].sample(
            n=len(positive_train_df))
        ])

    model = Class2NeuralNetwork(
        random_negative_rate=0.2,
        layer_sizes=[4],
        allele_positionwise_embedding_size=4,
        patience=5,
        peptide_convolutions=[
            {'kernel_size': 3, 'filters': 2, 'activation': "relu"},
        ],
        peptide_encoding={
            'vector_encoding_name': 'BLOSUM62',
            'alignment_method': 'right_pad',
            'max_length': 20,
        },
    )

    yield from train_and_check(train_df, model, alpha_sequences, beta_sequences)


def test_combination():
    # Can we generalize to an unseen allele?

    # Fake pseudosequences
    alpha_sequences = {
        "HLA-DRA*01:01": "AAAN",
    }
    beta_sequences = {
        "HLA-DRB1*01:01": "AAAA",
        "HLA-DRB1*03:01": "CAAA",
        "HLA-DRB1*04:01": "AAAC",
        "HLA-DRB1*05:01": "CAAC",
    }
    motifs = {
        "HLA-DRB1*01:01": "K.KK",
        "HLA-DRB1*03:01": "Q.KK",
        "HLA-DRB1*04:01": "K.KQ",
        "HLA-DRB1*05:01": "Q.KQ",
    }

    df = pandas.DataFrame(
        {"peptide": random_peptides(500000, length=15)}
    ).set_index("peptide")

    for (allele, motif) in motifs.items():
        df[allele] = (df.index.str.contains(motif)).astype(int)

    # Resample to have 1:1 binder / non-binder
    positive_train_df = df.loc[df.max(1) > 0.8]
    df = pandas.concat([
        positive_train_df,
        df.loc[~df.index.isin(positive_train_df.index)].sample(
            n=int(len(positive_train_df) / df.shape[1]))
    ])

    model = Class2NeuralNetwork(
        random_negative_rate=0.0,
        layer_sizes=[4],
        allele_positionwise_embedding_size=4,
        patience=10,
        peptide_convolutions=[
            {'kernel_size': 4, 'filters': 4, 'activation': "relu"},
        ],
        max_epochs=300,
        peptide_encoding={
            'vector_encoding_name': 'BLOSUM62',
            'alignment_method': 'right_pad',
            'max_length': 15,
        },
    )

    train_df = df.sample(frac=0.8).copy()
    train_df["HLA-DRB1*05:01"] = numpy.nan

    yield from train_and_check(
        df, model, alpha_sequences, beta_sequences, train_df=train_df)


def train_and_check(df, model, alpha_sequences, beta_sequences, train_df=None):
    print("Binders")
    print((df > 0.8).sum())

    print("Binder rate")
    print((df > 0.8).mean())

    if train_df is None:
        train_df = df.sample(frac=0.5)
    test_df = df.loc[~df.index.isin(train_df.index)]

    stacked = train_df.stack().reset_index().dropna()
    stacked.columns = ['peptide', 'allele', 'measurement_value']

    allele_encoding = make_allele_encoding_pair(
        stacked.allele, alpha_sequences, beta_sequences)

    print(model.hyperparameters)

    model.fit(
        stacked.peptide.values,
        affinities=stacked["measurement_value"].values,
        allele_encoding_pair=allele_encoding
    )

    yield from check_accuracy(
        train_df, model, alpha_sequences, beta_sequences, message="TRAIN")
    yield from check_accuracy(
        test_df, model, alpha_sequences, beta_sequences, message="TEST")


def check_accuracy(df, network, alpha_sequences, beta_sequences, message=""):
    stacked = df.stack().reset_index().dropna()
    stacked.columns = ['peptide', 'allele', 'measurement_value']

    allele_encoding = make_allele_encoding_pair(
        stacked.allele, alpha_sequences, beta_sequences)
    stacked["prediction"] = network.predict(
        stacked.peptide, allele_encoding_pair=allele_encoding)

    # Overall AUC
    stacked["binder"] = stacked.measurement_value > 0.8
    auc = roc_auc_score(stacked.binder, stacked.prediction)
    print(message, "Overall AUC", auc)
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
    print(message, "allele specific AUCs:")
    print(allele_specific_aucs)


    #import ipdb ; ipdb.set_trace()

