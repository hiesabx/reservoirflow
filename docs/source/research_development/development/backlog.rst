Backlog 📋
==========

Here we list the most important tasks that we need to work on. You may support us by working on one of these tasks below. You do not need to apply for a request; you can immediately fork our repository and start working on a new feature of your interest. As soon as you finish your work, then please make a pull request, so we can look at your work. Remember that you need to stick with the `Code of Conduct <../contribution.html#code-of-conduct>`_ and provide enough documentation and testing with your additional feature.

New Modules
-----------

.. attention:: 
  These modules have not been added yet. As soon as we have a added something 
  related, we will announce that in our 
  `Release Notes </release_notes/release_notes.html>`_ and 
  `News </community/news/news.html>`_.

- ``reservoirflow.pvt``: a module for PVT functionality. 
- ``reservoirflow.eos``: a module for EOS functionality. 
- ``reservoirflow.kr``: a module for relative permeability functionality.
- ``reservoirflow.pta``: a module for pressure transient analysis functionality. 
- ``reservoirflow.rta``: a module for rate transient analysis functionality.

Improving Current Modules
-------------------------

``reservoirflow.fluids``
^^^^^^^^^^^^^^^^^^^^^^^^
- `TwoPhase </api/reservoirflow.fluids.TwoPhase.html>`_ Fluid.
- `ThreePhase </api/reservoirflow.fluids.ThreePhase.html>`_ Fluid.
- `MultiPhase </api/reservoirflow.fluids.MultiPhase.html>`_ Fluid.

``reservoirflow.grids``
^^^^^^^^^^^^^^^^^^^^^^^
- `Radial </api/reservoirflow.grids.Radial.html>`_ Grid.
- `IrregularCartesian </api/reservoirflow.grids.IrregularCartesian.html>`_ Grid.

``reservoirflow.models``
^^^^^^^^^^^^^^^^^^^^^^^^
- `CompositionalModel </api/reservoirflow.models.CompositionalModel.html>`_ Model.
- `ThermalModel </api/reservoirflow.models.ThermalModel.html>`_ Model.

``reservoirflow.solutions``
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- `analytical </api/reservoirflow.solutions.analytical.html>`_ Solutions:
  `D1P1Solution  </api/reservoirflow.solutions.analytical.D1P1Solution.html>`_, 
  `D1P2Solution  </api/reservoirflow.solutions.analytical.D1P2Solution.html>`_, 
  `D1P3Solution  </api/reservoirflow.solutions.analytical.D1P3Solution.html>`_, 
  `D2P1Solution  </api/reservoirflow.solutions.analytical.D2P1Solution.html>`_, 
  `D2P2Solution  </api/reservoirflow.solutions.analytical.D2P2Solution.html>`_, 
  `D2P3Solution  </api/reservoirflow.solutions.analytical.D2P3Solution.html>`_, 
  `D3P1Solution  </api/reservoirflow.solutions.analytical.D3P1Solution.html>`_, 
  `D3P2Solution  </api/reservoirflow.solutions.analytical.D3P2Solution.html>`_, 
  `D3P3Solution  </api/reservoirflow.solutions.analytical.D3P3Solution.html>`_.
- `neurical </api/reservoirflow.solutions.neurical.html>`_ Solutions:
  `PinnSolution </api/reservoirflow.solutions.neurical.PinnSolution.html>`_, 
  `DeepONetSolution </api/reservoirflow.solutions.neurical.DeepONetSolution.html>`_.
- `numerical </api/reservoirflow.solutions.numerical.html>`_ Solutions:
  `FvmSolution </api/reservoirflow.solutions.numerical.FvmSolution.html>`_, 
  `FemSolution </api/reservoirflow.solutions.numerical.FemSolution.html>`_.

``reservoirflow.wells``
^^^^^^^^^^^^^^^^^^^^^^^
- `MultiCellWell </api/reservoirflow.wells.MultiCellWell.html>`_ Well.
- `DirectionalWell </api/reservoirflow.wells.DirectionalWell.html>`_ Well.

``reservoirflow.backends``
^^^^^^^^^^^^^^^^^^^^^^^^^^
- `Backend </api/reservoirflow.backends.Backend.html>`_.
- `NumPy </api/reservoirflow.backends.NumPy.html>`_.
- `PyTorch </api/reservoirflow.backends.PyTorch.html>`_.
- `TensorFlow </api/reservoirflow.backends.TensorFlow.html>`_.
- `JAX </api/reservoirflow.backends.JAX.html>`_.

``reservoirflow.utils``
^^^^^^^^^^^^^^^^^^^^^^^
- Import/Export Eclipse.
- Import/Export tNavigator.
- Unit converter.

|

.. include:: /_static/comments_section.rst