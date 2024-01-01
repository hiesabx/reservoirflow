from reservoirflow.backends import Backend


class JAX(Backend):
    """JAX backend class.

    .. caution::
        This class is not available.

    Returns
    -------
    Backend
        Backend object.

    References
    ----------
    `JAX <https://jax.readthedocs.io/en/latest/index.html>`_.
    """

    name = "JAX"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")

    def transpose(self):
        raise NotImplementedError

    def ones(self):
        raise NotImplementedError
