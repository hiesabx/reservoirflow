# ReservoirFlow: a Petroleum Reservoir Simulation Library in Python

![five_spot_single_phase](/docs\source\user_guide\tutorials\tutorial_five_spot_single_phase\grid_animated.gif)\
*Example: Pressure Distribution of Single Phase Flow in Five Spot Wells Patterns.*

**Table of Content:**

- [Introduction](#introduction)
- [Installation](#installation)
- [Import Convention](#import-convention)
- [Version](#version)
- [License](#license)
- [Contact](#contact)

## Introduction

*ReservoirFlow* <`reservoirflow`> is a modern open-source Python library developed by [Bayanatics](https://github.com/zakgrin) to study and model the phenomenon of fluid flow in porous media (i.e. Reservoir Simulation) using different approaches including analytical solutions, numerical solutions, and machine learning solutions. *ReservoirFlow* is designed based on the modern Python stack for data science, scientific computing, and machine learning with the objective to support high-performance computing including multithreading, parallelism, GPU, and TPU.

*ReservoirFlow* brings reservoir simulation to the Python ecosystem including data analytics and machine learning tools to empower intelligent fields where AI specialists can deploy their models in containers that will be ready to make real-time optimization for any well in the field. Having this option allows for better integration and support with intelligent fields. In contrast to commercial black-box software where reservoir simulation studies are relatively isolated, important actions can be immediately predicted and made available for the field hardware to execute.

*ReservoirFlow* aims to achieving a high quality open research and science for reservoir engineering and simulation. Solutions that can combine the strength of scientific computing with the power of machine learning (e.g. deep learning) for different objectives (e.g. reverse computing, interpolation or extrapolation, etc.). Below are few examples of the problems that will be tackled in the future:

- Real-time reservoir management and production optimization using Cloud Computing and IoT.
- Reinforcement learning to achieve better production strategies for specific goals (e.g. maximize recovery, accelerate production).
- History matching using machine learning.
- Advanced computing such as GPU, TPU and Quantum Computing.
- Scientific Machine learning using Physics-informed neural networks (PINNs) or DeepONets.

An open-source reservoir simulation library within the Python ecosystem is also very important to students, universities, researchers, engineers, and practitioners. Therefore, the growth of this tool can only be taken as a positive growth for a new and healthy oil & gas community that we try to create. However, this requires a huge support to meet the upcoming challenges that we are looking for, see [Support Us](/support_us.html).

## Installation

Install `reservoirflow` directly from [PyPi](https://pypi.org/):

```bash
pip install reservoirflow
```

For more information about the installation process, see: [Getting Started](/user_guide/getting_started/getting_started.html) in the documentation.

## Import Convention

The following convention is used to import `reservoirflow` after installation:

```python
import reservoirflow as rf
```

The abbreviation `rf` refers to `reservoirflow` where all modules under this library can be accessed. `rf` is also used throughout the documentation. We recommend our users to stick with this convention.

## Version

[Semantic Versioning](https://semver.org/) is used for the version numbers. Since this library is still under active development, `major=0` is used until the first stable version is reached. The first version `v0.1.0` was released on January 1, 2024. The current version is `v0.1.0`. To know about which features are currently supported, check [Capabilities](capabilities.html).

**Version History:**

||**Version**|**Status**|**Release Date (dd.mm.yyyy)**|
|-|-|-|-|
|1|`v0.1.0`|current version|01.01.2024|
|2|`v0.1.1`|*under development*|*ongoing*|

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). Detailed license can be found [here](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode). The current license does not allow commercial use. For commercial applications, please [contact](#contact) us.

## Contact

You can contact us directly on the following email: <reservoirflow.contact@gmail.com>.
