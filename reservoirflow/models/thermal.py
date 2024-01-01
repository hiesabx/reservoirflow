from reservoirflow.models.model import Model


class Thermal(Model):
    """Thermal model class.

    .. caution::
        This class is not available.

    Returns
    -------
    Model
        Model object.
    """

    name = "Thermal"

    def __init__(self, **kwargs):
        raise NotImplementedError("This class is not implemented.")
