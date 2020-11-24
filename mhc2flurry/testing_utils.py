"""
Utilities used in MHCflurry unit tests.
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
    K.clear_session()
