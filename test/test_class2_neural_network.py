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

from mhc2flurry.class2_neural_network import Class2NeuralNetwork
from mhc2flurry.downloads import get_path
from mhc2flurry.common import random_peptides

from mhc2flurry.testing_utils import cleanup, startup
teardown = cleanup
setup = startup


def test_memorize():
    train_df = pandas.DataFrame({"peptide": random_peptides(1000, length=15)})
    train_df["HLA-DRB*01:01"] = 20000
    train_df["HLA-DRB*03:01"] = 20000
    train_df.loc[
        train_df.peptide.str.contains("AAK"), "HLA-DRB*01:01"
    ] = 50.0
    train_df.loc[
        train_df.peptide.str.contains("QAQ"), "HLA-DRB*03:01"
    ] = 50.0

    # Memorize the dataset.
    model = Class2NeuralNetwork(
        random_negative_rate=1.0,
        layer_sizes=[8],
        patience=5,
        peptide_convolutions=[
            {'kernel_size': 3, 'filters': 8, 'activation': "relu"},
            {'kernel_size': 1, 'filters': 4, 'activation': "relu"},
            {'kernel_size': 4, 'filters': 4, 'activation': "relu"},
        ],
    )
    print(model.hyperparameters)

    import ipdb ; ipdb.set_trace()

    model.fit(
        use_train_df.peptide.values,
        affinities=use_train_df["measurement_value"].values,
        inequalities=use_train_df["measurement_inequality"].values,
        allele_encoding_pair=allele_encoding_pair
    )
