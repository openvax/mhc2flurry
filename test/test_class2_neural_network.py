import logging
logging.getLogger('tensorflow').disabled = True
logging.getLogger('matplotlib').disabled = True

import numpy
import tensorflow.random
numpy.random.seed(0)
tensorflow.random.set_seed(0)

import pandas
from sklearn.metrics import roc_auc_score

import mhcgnomes

from mhc2flurry.allele_encoding_pair import AlleleEncodingPair
from mhc2flurry.allele_encoding import AlleleEncoding
from mhc2flurry.class2_neural_network import Class2NeuralNetwork
from mhc2flurry.common import random_peptides

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


def test_simple():
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
        minibatch_size=1024,
        random_negative_rate=1.0,
        layer_sizes=[4],
        allele_positionwise_embedding_size=4,
        patience=10,
        max_epochs=500,
        peptide_convolutions=[
            {'kernel_size': 3, 'filters': 8, 'activation': "relu"},
        ],
        peptide_encoding={
            'vector_encoding_name': 'BLOSUM62',
            'alignment_method': 'right_pad',
            'max_length': 20,
        },
    )

    train_and_check(train_df, model, alpha_sequences, beta_sequences)


def test_combination():
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
        "HLA-DRB1*01:01": "K.AK",
        "HLA-DRB1*03:01": "Q.CK",
        "HLA-DRB1*04:01": "K.DQ",
        "HLA-DRB1*05:01": "Q.EQ",
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
        minibatch_size=1024,
        random_negative_rate=1.0,
        layer_sizes=[4],
        allele_positionwise_embedding_size=4,
        patience=10,
        peptide_convolutions=[
            {'kernel_size': 4, 'filters': 12, 'activation': "relu"},
        ],
        max_epochs=500,
        peptide_encoding={
            'vector_encoding_name': 'BLOSUM62',
            'alignment_method': 'right_pad',
            'max_length': 15,
        },
    )

    train_df = df.sample(frac=0.8).copy()

    # Can we generalize to an unseen allele?
    # So far, haven't gotten this to work, so leaving this line commented.
    #train_df["HLA-DRB1*05:01"] = numpy.nan

    train_and_check(
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

    check_accuracy(
        train_df, model, alpha_sequences, beta_sequences, message="TRAIN")
    check_accuracy(
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
    assert auc > 0.7, message

    # Can we discern a binder for one allele from another?
    binder_peptides = stacked.loc[stacked.binder].peptide.unique()
    stacked_binders = stacked.loc[stacked.peptide.isin(binder_peptides)]
    allele_specific_aucs = []
    for (allele, sub_df) in stacked_binders.groupby("allele"):
        print(allele)
        print(sub_df)
        auc = roc_auc_score(sub_df.binder.values, sub_df.prediction.values)
        allele_specific_aucs.append((allele, auc))

    allele_specific_aucs = pandas.DataFrame(
        allele_specific_aucs, columns=["allele", "auc"])
    print(message, "allele specific AUCs:")
    print(allele_specific_aucs)

    print(message, "Mean predictions")
    print(stacked_binders.groupby(["allele", "binder"]).prediction.mean())

    for _, row in allele_specific_aucs.iterrows():
        assert row.auc > 0.8, (message, row.allele)


