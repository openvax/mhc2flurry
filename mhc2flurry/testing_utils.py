"""
Utilities used in MHC2flurry unit tests.
"""
from .common import configure_tensorflow


def startup():
    """
    Configure Keras backend for running unit tests.
    """
    configure_tensorflow("tensorflow-cpu", num_threads=2)


def cleanup():
    """
    Clear tensorflow session and other process-wide resources.
    """
    import tensorflow.keras.backend as K
    #Class2NeuralNetwork.clear_model_cache()
    K.clear_session()
