from reservoirflow.backends import Backend


class PyTorch(Backend):
    """PyTorch backend class.

    .. caution::
        This class is not available.

    Returns
    -------
    Backend
        Backend object.

    References
    ----------
    `PyTorch <https://pytorch.org/>`_.
    """

    name = "PyTorch"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def transpose(self):
        raise NotImplementedError

    def ones(self):
        raise NotImplementedError
