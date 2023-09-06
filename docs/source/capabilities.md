# Capabilities

A long-term plan with high ambitions is already in place for this library. As of now, this tool is still very small relative to the entire scope of work (i.e. only about 5%). Supported capabilities so far are shown below:

|**Feature**|**Type**|**Support**|**Starting From**|**Comment**|
|-----------|--------|-----------|-----------------|-----------|
|**Model Type**|Black Oil|YES|`v0.1.0`|Only single phase on Cartesian grid.|
||Compositional|No|`-`||
|**Flow Dimensions**|1D|Yes|`v0.1.0`|Flow dimension is defined based on grid shape.|
||2D|Yes|`v0.1.0`|Flow dimension is defined based on grid shape.|
||3D|Yes|`v0.1.0`|Flow dimension is defined based on grid shape.|
|**Grid Types**|Regular Cartesian|Yes|`v0.1.0`||
||Single Well Radial|No|`-`||
||Cartesian|No|`-`||
|**Fluid Phases**| Single Phase|Yes|`v0.1.0`||
||Two Phases|No|`-`||
||Three Phases|No|`-`||
|**Fluid Compressibility**| Incompressible |Yes|`v0.1.0`||
|| Slightly Compressible|Yes|`v0.1.0`||
|| Compressible|No|`-`||
|**Rock Compressibility**| Incompressible |Yes|`v0.1.0`||
|| Slightly Compressible|Yes|`v0.1.0`||
|**Experiments**|Core-Flooding|No|`-`||
|| Slim-Tube|No|`-`||
|**Well Types**|Single Cell|Yes|`v0.1.0`||
||Multiple Cells|No|`-`||
||Horizontal|No|`-`||
||Slanted|No|`-`||
|**History Matching**|Conventional|No| `-`||
||Machine Learning|No|`-`||
|**Production Optimization**|Reinforcement Learning|No|`-`||
|**Solutions**|Numerical [FDM]|Yes|`v0.1.0`|Requires an `Initializer` and a `Solver`|
||Analytical [1D]|Yes|`v0.1.0`|Available only for `1D Single Phase` problems.|
||Neurical [1D]|Yes|`v0.1.0`|Neural networks based solution (e.g., PINNs).|
|**Numerical Initializers**|Vectorized|Yes|`v0.1.0`|Used for efficient computing.|
||Symbolic|Yes|`v0.1.0`|Used for debugging.|
|**Numerical Solvers**|Iterative Solvers|Yes|`v0.1.0`|Requires sparse matrices.|
||Neurical Solvers*|No|`-`|Solving a system of linear equations using neural networks.|
|**Efficient Computing**|Sparse Matrices|Yes|`v0.1.0`|Required for iterative solvers.|
||Threading|Yes|`v0.1.0`|Used for numerical symbolic initializers.|
||Concurrent|Yes|`v0.1.0`|Used for numerical symbolic initializers.|
||GPU Computing|No|`-`||
||TPU Computing|No|`-`||
||Quantum Computing|No|`-`||
|

(*) Innovative solution introduced within this work.
