# Capabilities

A long-term plan with high ambitions is already in place for this library. As of now, this tool is still very small relative to the entire scope of work (i.e. only about 5%). Supported capabilities so far are shown below:

|**Feature**                |**Type**               |**Support**|**Starting From**  |**Comment**|
| -------------------------:| --------------------- | --------- | ----------------- | --------- |
|**Model Type**             |Black Oil              |yes        |`v0.1.0`           |Only single phase on regular cartesian grid.|
|                           |Compositional          |no         |`-`                ||
|                           |Thermal                |no         |`-`                ||
|**Flow Dimensions**        |1D                     |yes        |`v0.1.0`           |Flow dimension is defined based on grid shape.|
|                           |2D                     |yes        |`v0.1.0`           |Flow dimension is defined based on grid shape.|
|                           |3D                     |yes        |`v0.1.0`           |Flow dimension is defined based on grid shape.|
|**Grid Types**             |Regular Cartesian      |yes        |`v0.1.0`           ||
|                           |Radial                 |no         |`-`                ||
|                           |Irregular Cartesian    |no         |`-`                ||
|**Fluid Phases**           |Single Phase           |yes        |`v0.1.0`           ||
|                           |Two Phases             |no         |`-`                ||
|                           |Three Phases           |no         |`-`                ||
|                           |Compositional          |no         |`-`                ||
|**Fluid Compressibility**  |Incompressible         |yes        |`v0.1.0`           ||
|                           |Slightly Compressible  |yes        |`v0.1.0`           ||
|                           |Compressible           |no         |`-`                ||
|**Rock Compressibility**   |Incompressible         |yes        |`v0.1.0`           ||
|                           |Slightly Compressible  |yes        |`v0.1.0`           ||
|**Experiments**            |Core-Flooding          |no         |`-`                ||
|                           |Slim-Tube              |no         |`-`                ||
|**Well Types**             |Single Cell            |yes        |`v0.1.0`           ||
|                           |Multiple Cells         |no         |`-`                ||
|                           |Directional            |no         |`-`                ||
|**History Matching**       |Conventional           |no         |`-`                ||
|                           |Machine Learning       |no         |`-`                ||
|**Production Optimization**|Reinforcement Learning |no         |`-`                ||
|**Solutions**              |Numerical [FDM]        |yes        |`v0.1.0`           |Requires an `Initializer` and a `Solver`|
|                           |Analytical [1D]        |yes        |`v0.1.0`           |Available only for `1D Single Phase` problems.|
|                           |Neurical [1D]          |yes        |`v0.1.0`           |Neural networks based solution (e.g., PINNs).|
|**Numerical Initializers** |Vectorized             |yes        |`v0.1.0`           |Used for efficient computing.|
|                           |Symbolic               |yes        |`v0.1.0`           |Used for debugging.|
|**Numerical Solvers**      |Iterative Solvers      |yes        |`v0.1.0`           |Requires sparse matrices.|
|                           |Neurical Solvers*      |no         |`-`                |Solving a system of linear equations using neural networks.|
|**Efficient Computing**    |Sparse Matrices        |yes        |`v0.1.0`           |Required for iterative solvers.|
|                           |Threading              |yes        |`v0.1.0`           |Used for numerical symbolic initializers.|
|                           |Concurrent             |yes        |`v0.1.0`           |Used for numerical symbolic initializers.|
|                           |GPU Computing          |no         |`-`                ||
|                           |TPU Computing          |no         |`-`                ||
|                           |Quantum Computing      |no         |`-`                ||
|

(*) Innovative solution introduced within this work.

Feel free to make a comment below.

```{include} /_static/comments_section.md
```
