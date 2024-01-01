from reservoirflow.backends import Backend


class NumPy(Backend):
    """NumPy backend class.

    .. caution::
        This class is not available.

    Returns
    -------
    Backend
        Backend object.

    References
    ----------
    `NumPy <https://numpy.org/doc/stable/index.html>`_.
    """

    name = "NumPy"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def transpose(self):
        raise NotImplementedError

    def ones(self):
        raise NotImplementedError
