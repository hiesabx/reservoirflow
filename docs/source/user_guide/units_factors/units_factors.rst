Units & Factors 📐
==================

ReservoirFlow has a fixed ``unit`` system used when an object is initiated (i.e., object for the modules such as: `grids </api/reservoirflow.grids.html>`_, `fluids </api/reservoirflow.fluids.html>`_, `models </api/reservoirflow.models.html>`_). The value of ``unit`` property can be set to ``"field"``, ``"metric"``, or ``"lab"``. By default, ``unit="field"`` is used. This system was inspired by :footcite:t:`abou2013petroleum`. This system was inspired by :cite:t:`aziz2002petroleum` and :cite:t:`abou2013petroleum`.

.. warning::
   As of now, changing class unit will not automatically convert values. In the near future, a unit converter will be developed to match units and easily convert units from one system to another.


Units
-----

.. include:: units_table.rst

Factors
-------

.. include:: factors_table.rst

|

.. toctree::
   :maxdepth: 1
   :caption: See Also

   access_units_factors

|

.. include:: /_static/references_section.rst

|

.. include:: /_static/comments_section.rst
