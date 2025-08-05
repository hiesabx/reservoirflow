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
- `Compositional </api/reservoirflow.models.Compositional.html>`_ Model.
- `Thermal </api/reservoirflow.models.Thermal.html>`_ Model.

``reservoirflow.solutions``
^^^^^^^^^^^^^^^^^^^^^^^^^^^
- `analytical </api/reservoirflow.solutions.analytical.html>`_ Solutions:
  `D1P1  </api/reservoirflow.solutions.analytical.D1P1.html>`_, 
  `D1P2  </api/reservoirflow.solutions.analytical.D1P2.html>`_, 
  `D1P3  </api/reservoirflow.solutions.analytical.D1P3.html>`_, 
  `D2P1  </api/reservoirflow.solutions.analytical.D2P1.html>`_, 
  `D2P2  </api/reservoirflow.solutions.analytical.D2P2.html>`_, 
  `D2P3  </api/reservoirflow.solutions.analytical.D2P3.html>`_, 
  `D3P1  </api/reservoirflow.solutions.analytical.D3P1.html>`_, 
  `D3P2  </api/reservoirflow.solutions.analytical.D3P2.html>`_, 
  `D3P3  </api/reservoirflow.solutions.analytical.D3P3.html>`_.
- `neurical </api/reservoirflow.solutions.neurical.html>`_ Solutions:
  `PINN </api/reservoirflow.solutions.neurical.PINN.html>`_, 
  `DeepONet </api/reservoirflow.solutions.neurical.DeepONet.html>`_.
- `numerical </api/reservoirflow.solutions.numerical.html>`_ Solutions:
  `FVM </api/reservoirflow.solutions.numerical.FVM.html>`_, 
  `FEM </api/reservoirflow.solutions.numerical.FEM.html>`_.

``reservoirflow.wells``
^^^^^^^^^^^^^^^^^^^^^^^
- `MultiCell </api/reservoirflow.wells.MultiCell.html>`_ Well.
- `Directional </api/reservoirflow.wells.Directional.html>`_ Well.

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