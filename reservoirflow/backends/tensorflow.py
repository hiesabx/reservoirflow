from reservoirflow.backends import Backend


class TensorFlow(Backend):
    """TensorFlow backend class.

    .. caution::
        This class is not available.

    Returns
    -------
    Backend
        Backend object.

    References
    ----------
    `TensorFlow <https://www.tensorflow.org/>`_.
    """

    name = "TensorFlow"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def transpose(self):
        raise NotImplementedError

    def ones(self):
        raise NotImplementedError
