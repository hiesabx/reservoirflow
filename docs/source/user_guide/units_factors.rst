***************
Units & Factors
***************
Version: |release|

ReservoirFlow has a fixed ``unit`` system used when a class is initiated (e.g., Grid: `RegularCartesian </api/reservoirflow.grids.html#reservoirflow.grids.RegularCartesian>`_, Fluid: `SinglePhase </api/reservoirflow.fluids.html#reservoirflow.fluids.SinglePhase>`_, Model: `BlackOil </api/reservoirflow.models.html#reservoirflow.models.BlackOil>`_). The value of ``unit`` property can be set to ``"field"``, ``"metric"``, or ``"lab"``. By default, ``unit="field"`` is used. As of now, changing class unit will not automatically convert values. In the near future, a unit converter will be developed to match units and easily convert units from one type to another.

Units
#####

.. include:: UNITS.rst

Factors
#######

.. include:: FACTORS.rst