# Capabilities

A long-term plan with high ambitions is already in place for this library. As of now, this tool is still very small relative to the entire scope of work (i.e. only about 5%). Planned and supported capabilities are shown below:

|**Module**                |**SubModule/Class**      |**Support**|**Starting From**  |**Comment**|
|:-------------------------:|------------------------|-----------|-------------------|-----------|
|**[models]**               |[BlackOil]              |yes        |`v0.1.0`           |Only single phase on regular cartesian grid.|
|                           |[Compositional]         |no         |`-`                ||
|                           |[Thermal]               |no         |`-`                ||
|**[grids]**                |[RegularCartesian]      |yes        |`v0.1.0`           |Includes: rock compressibility.|
|                           |[IrregularCartesian]    |no         |`-`                ||
|                           |[Radial]                |no         |`-`                ||
|**[fluids]**               |[SinglePhase]           |yes        |`v0.1.0`           |Includes: fluid compressibility.|
|                           |[TwoPhase]              |no         |`-`                ||
|                           |[ThreePhase]            |no         |`-`                ||
|                           |[MultiPhase]            |no         |`-`                ||
|**[wells]**                |[SingleCell]            |yes        |`v0.1.0`           ||
|                           |[MultiCell]             |no         |`-`                ||
|                           |[Directional]           |no         |`-`                ||
|**[scalers]**              |[Dummy]                 |yes        |`v0.1.0`           ||
|                           |[MinMax]                |yes        |`v0.1.0`           ||
|**[solutions]**            |[numerical]             |yes        |`v0.1.0`           |Includes: Finit-Difference-Method (`FDM`)|
|                           |[analytical]            |yes        |`v0.1.0`           |Includes: 1-Dimentional-1Phase (`D1P1`)|
|                           |[neurical]              |yes        |`v0.1.0`           |Includes: Physics-Informed-Neural-Network (`PINN`)|
|**Experiments**            |Core-Flooding           |no         |`-`                ||
|                           |Slim-Tube               |no         |`-`                ||
|**History Matching**       |Conventional            |no         |`-`                ||
|                           |Machine Learning        |no         |`-`                ||
|**Production Optimization**|Reinforcement Learning  |no         |`-`                ||
|**Numerical Initializers** |Vectorized              |yes        |`v0.1.0`           |Used for efficient computing.|
|                           |Symbolic                |yes        |`v0.1.0`           |Used for debugging.|
|**Numerical Solvers**      |Iterative Solvers       |yes        |`v0.1.0`           |Requires sparse matrices.|
|                           |Neurical Solvers*       |no         |`-`                |Solving a system of linear equations using neural networks.|
|**Efficient Computing**    |Sparse Matrices         |yes        |`v0.1.0`           |Required for iterative solvers.|
|                           |Threading               |yes        |`v0.1.0`           |Used for numerical symbolic initializers.|
|                           |Concurrent              |yes        |`v0.1.0`           |Used for numerical symbolic initializers.|
|                           |GPU Computing           |no         |`-`                ||
|                           |TPU Computing           |no         |`-`                ||
|                           |Quantum Computing       |no         |`-`                ||
 
(*) Innovative solution introduced within this work.


[models]: /api/reservoirflow.models.html
[BlackOil]: /api/reservoirflow.models.BlackOil.html
[Compositional]: /api/reservoirflow.models.Compositional.html
[Thermal]: /api/reservoirflow.models.Thermal.html

[grids]: /api/reservoirflow.grids.html
[RegularCartesian]: /api/reservoirflow.grids.RegularCartesian.html
[Radial]: /api/reservoirflow.grids.Radial.html
[IrregularCartesian]: /api/reservoirflow.grids.IrregularCartesian.html

[fluids]: /api/reservoirflow.fluids.html
[SinglePhase]: /api/reservoirflow.fluids.SinglePhase.html
[TwoPhase]: /api/reservoirflow.fluids.TwoPhase.html
[ThreePhase]: /api/reservoirflow.fluids.ThreePhase.html
[MultiPhase]: /api/reservoirflow.fluids.MultiPhase.html

[wells]: /api/reservoirflow.wells.html
[SingleCell]: /api/reservoirflow.wells.SingleCell.html
[MultiCell]: /api/reservoirflow.wells.MultiCell.html
[Directional]: /api/reservoirflow.wells.Directional.html

[scalers]: /api/reservoirflow.scalers.html
[Dummy]: /api/reservoirflow.scalers.SingleCell.html
[MinMax]: /api/reservoirflow.scalers.MultiCell.html

[solutions]: /api/reservoirflow.solutions.html
[analytical]: /api/reservoirflow.solutions.analytical.html
[neurical]: /api/reservoirflow.solutions.neurical.html
[numerical]: /api/reservoirflow.solutions.numerical.html


```{include} /_static/comments_section.md
```
